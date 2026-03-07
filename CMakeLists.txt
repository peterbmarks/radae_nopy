cmake_minimum_required(VERSION 3.16)
project(radae_nopy C)

set(CMAKE_C_STANDARD 11)

# Build opus with FARGAN/LPCNet support
option(AVX "Enable AVX CPU optimizations." ON)
include(cmake/BuildOpus.cmake)

add_subdirectory(src)

# C Port of Core Encoder/decoder ------------------------------------------------------------------------------------

if (NOT WIN32)
# we test by comparing loss of features_in/features_out, can all happen at rate Rs.  We load model05.bin weights as 
# compiled-in weights are for model19_check3
add_test(NAME c_encoder_model5
         COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; 
         ${CMAKE_CURRENT_BINARY_DIR}/src/lpcnet_demo -features wav/brian_g8sez.wav features_in.f32; \
         cat features_in.f32 | ${CMAKE_CURRENT_BINARY_DIR}/src/test_rade_enc 1 0 ${CMAKE_SOURCE_DIR}/bin/model05.bin > z_c.f32; \
         python3 stateful_encoder.py model05/checkpoints/checkpoint_epoch_100.pth features_in.f32 /dev/null --read_latent z_c.f32 --loss_test 0.2")
         set_tests_properties(c_encoder_model5 PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# we need the DSP code in the loop to test the core encoder with model19_check3 as the bottleneck (3) is at rate Fs
add_test(NAME c_encoder_model19_check3
         COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; 
         ${CMAKE_CURRENT_BINARY_DIR}/src/lpcnet_demo -features wav/brian_g8sez.wav features_in.f32; \
         cat features_in.f32 | python3 radae_txe.py | python3 radae_rxe.py > features_out.f32; \
         cat features_in.f32 | ${CMAKE_CURRENT_BINARY_DIR}/src/test_rade_enc 3 1 > z_c.f32; \
         cat z_c.f32 | python3 radae_txe.py --bypass_enc | python3 radae_rxe.py > features_rx_out.f32; \
         python3 loss.py features_in.f32 features_out.f32 --features_hat2 features_rx_out.f32 --compare")
         set_tests_properties(c_encoder_model19_check3 PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# compare Python decoder (features_out.f32) C decoder (features_c.f32), for same latent inputs z.f32
# note inference.sh creates features_in.f32 & features_out.f32
add_test(NAME c_decoder_model5
         COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; 
         ./inference.sh model05/checkpoints/checkpoint_epoch_100.pth wav/brian_g8sez.wav /dev/null --write_latent z.f32; \
         cat z.f32 | ${CMAKE_CURRENT_BINARY_DIR}/src/test_rade_dec 0 ${CMAKE_SOURCE_DIR}/bin/model05.bin > features_c.f32; \
         python3 loss.py features_in.f32 features_out.f32 --features_hat2 features_c.f32 --compare")
         set_tests_properties(c_decoder_model5 PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# As well as testing C core decoder, this test also uses radae_rxe.py in bypass mode
add_test(NAME c_decoder_model19_check3
         COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; 
                       ./inference.sh model19_check3/checkpoints/checkpoint_epoch_100.pth wav/brian_g8sez.wav /dev/null \
                       --rate_Fs --pilots --pilot_eq --eq_ls --cp 0.004 --bottleneck 3 --auxdata --correct_freq_offset \
                       --write_rx rx.f32; \
         cat rx.f32 | python3 radae_rxe.py --bypass_dec > z_hat.f32
         cat z_hat.f32 | ${CMAKE_CURRENT_BINARY_DIR}/src/test_rade_dec 1 > features_c.f32; \
         python3 loss.py features_in.f32 features_out.f32 --features_hat2 features_c.f32 --compare")
         set_tests_properties(c_decoder_model19_check3 PROPERTIES PASS_REGULAR_EXPRESSION "PASS")
endif(NOT WIN32)

# Test embedded data (--auxdata) use for false sync detection, with C core decoder, with C callable embedded Python API.
# --foff_err forces a false sync state after first sync. See also "radae_rx_aux_mpp" above.  Tests uw_errors are being passed
# back from C decoder to Python state machine
add_test(NAME c_decoder_aux_mpp
        COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; \
                       test/make_g.sh; \
                       ./inference.sh model19_check3/checkpoints/checkpoint_epoch_100.pth wav/all.wav /dev/null \
                       --rate_Fs --pilots --pilot_eq --eq_ls --cp 0.004 --bottleneck 3 --time_offset -16 --auxdata \
                       --EbNodB 4 --freq_offset -11 --g_file g_mpp.f32 --write_rx rx.f32 \
                       --prepend_noise 1 --append_noise 3 --end_of_over --correct_freq_offset; \
                       cat rx.f32 | PYTHONPATH='.' ${CMAKE_CURRENT_BINARY_DIR}/src/radae_rx 1 > features_c.f32; \
                       python3 loss.py features_in.f32 features_c.f32 --loss 0.3 --clip_start 300")
                       set_tests_properties(c_decoder_aux_mpp PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# EOO data -------------------------------------------------------------------------------------------
# Pythion Tx & Rx, no noise
add_test(NAME radae_eoo_data_py
        COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; \
                       ${CMAKE_CURRENT_BINARY_DIR}/src/lpcnet_demo -features wav/brian_g8sez.wav features_in.f32; \
                       cat features_in.f32 | python3 radae_txe.py --eoo_data_test > rx.f32; \
                       cat rx.f32 | python3 radae_rxe.py -v 2 --eoo_data_test > /dev/null")
                       set_tests_properties(radae_eoo_data_py PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# C tx & rx, no noise. Note Python radae_txe.py just generates eoo_tx.f32 for C radae_tx
add_test(NAME radae_eoo_data_c
        COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; \
                       ${CMAKE_CURRENT_BINARY_DIR}/src/lpcnet_demo -features wav/brian_g8sez.wav features_in.f32; \
                       cat features_in.f32 | python3 radae_txe.py --eoo_data_test > /dev/null; \
                       cat features_in.f32 | PYTHONPATH='.' ${CMAKE_CURRENT_BINARY_DIR}/src/radae_tx > rx.f32; \
                       cat rx.f32 | PYTHONPATH='.' ${CMAKE_CURRENT_BINARY_DIR}/src/radae_rx > /dev/null; \
                       python3 eoo_ber.py eoo_tx.f32 eoo_rx.f32")
                       set_tests_properties(radae_eoo_data_c PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# C tx & rx, over a multipath channel, we set up 5 "overs", each with an EOO chunk of data.  Just one of them needs
# to have a BER < 5% for a PASS.  About 6dB SNR
add_test(NAME radae_eoo_data_mpp
        COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; test/make_g.sh; \
                       ${CMAKE_CURRENT_BINARY_DIR}/src/lpcnet_demo -features wav/brian_g8sez.wav features_in.f32; \
                       cat features_in.f32 | python3 radae_txe.py --eoo_data_test > /dev/null; \
                       cat features_in.f32 | PYTHONPATH='.' ${CMAKE_CURRENT_BINARY_DIR}/src/radae_tx > tx.f32; \
                       cat tx.f32 tx.f32 tx.f32 tx.f32 tx.f32 > tx_2.f32;
                       cat tx_2.f32 | python3 f32toint16.py --real --scale 8192 > tx.raw; \
                       ${CODEC2_DEV_BUILD_DIR}/src/ch tx.raw rx.raw --No -24 --mpp --fading_dir .; \
                       cat rx.raw | python3 int16tof32.py --zeropad > rx.f32; \
                       cat rx.f32 | PYTHONPATH='.' ${CMAKE_CURRENT_BINARY_DIR}/src/radae_rx > /dev/null; \
                       python3 eoo_ber.py eoo_tx.f32 eoo_rx.f32")
                       set_tests_properties(radae_eoo_data_mpp PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# BBFM -----------------------------------------------------------------------------------------------

# single carrier modem internal (inside single_carrier class) tests
add_test(NAME bbfm_sc_internal
        COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; \
                       python3 -c 'from radae import single_carrier,single_carrier_tests; single_carrier_tests()'")
                       set_tests_properties(bbfm_sc_internal PROPERTIES PASS_REGULAR_EXPRESSION "ALL PASS")

# single carrier modem stand alone tx/rx, using BPSK symbol/BER test mode
add_test(NAME bbfm_sc_ber
        COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; \
                       ./bbfm_inference.sh model_bbfm_01/checkpoints/checkpoint_epoch_100.pth wav/brian_g8sez.wav /dev/null --write_latent z.f32; \
                       cat z.f32 | python3 sc_tx.py --ber_test > t.int16; \
                       cat t.int16 | python3 sc_rx.py --ber_test --target_ber 0 > /dev/null")
                       set_tests_properties(bbfm_sc_ber PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# single carrier modem stand alone tx/rx, using continuously valued BBFM symbols, compare loss in decoded features
add_test(NAME bbfm_sc_loss
        COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; \
                       ./bbfm_inference.sh model_bbfm_01/checkpoints/checkpoint_epoch_100.pth wav/brian_g8sez.wav /dev/null --write_latent z.f32; \
                       cat z.f32 | python3 sc_tx.py | python3 sc_rx.py > z_hat.f32; \
                       ./bbfm_rx.sh model_bbfm_01/checkpoints/checkpoint_epoch_100.pth z_hat.f32 /dev/null; \
                       python3 loss.py features_in.f32 features_out.f32 --features_hat2 features_rx_out.f32 --compare")
                       set_tests_properties(bbfm_sc_loss PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

# single carrier modem stand alone tx/rx with external 300-2700Hz band pass filter and no noise, measure loss.  In practice, SNR will be
# quite high, so channel distortions other than noise may dominate
add_test(NAME bbfm_sc_bpf_loss
        COMMAND sh -c "cd ${CMAKE_SOURCE_DIR}; \
                       ./bbfm_inference.sh model_bbfm_01/checkpoints/checkpoint_epoch_100.pth wav/brian_g8sez.wav /dev/null --write_latent z.f32; \
                       cat z.f32 | python3 sc_tx.py | ${CODEC2_DEV_BUILD_DIR}/src/ch - - | python3 sc_rx.py > z_hat.f32; \
                       ./bbfm_rx.sh model_bbfm_01/checkpoints/checkpoint_epoch_100.pth z_hat.f32 /dev/null; \
                       python3 loss.py features_in.f32 features_out.f32 --features_hat2 features_rx_out.f32 --compare")
                       set_tests_properties(bbfm_sc_bpf_loss PROPERTIES PASS_REGULAR_EXPRESSION "PASS")
