# Experimental version of FreeDV RADE without python

Based on work from [David Rowe](https://github.com/drowe67/radae) et al.

## Build
```
cd radae_nopy
mkdir build && cd build
cmake ..
make -j4
```
## Usage
### Encode: WAV → IQ
```
sox ../voice.wav -r 16000 -t .s16 -c 1 - | \
  ./src/lpcnet_demo -features /dev/stdin - | \
  ./src/radae_tx > tx.iq
```

### Decode: IQ → WAV  
```
cat tx.iq | \
  ./src/radae_rx | \
  ./src/lpcnet_demo -fargan-synthesis /dev/stdin - | \
  sox -t .s16 -r 16000 -c 1 - decoded.wav
```

## Files
| File                                          | Purpose                                             |
| --------------------------------------------- | --------------------------------------------------- |
| `src/rade_dsp.h/c`                            | Complex math utilities, constants, pilot generation |
| `src/rade_ofdm.h/c`                           | OFDM modulation/demodulation with DFT matrices      |
| `src/rade_bpf.h/c`                            | Complex bandpass filter                             |
| `src/rade_tx.h/c`                             | Transmitter (encoder + OFDM mod)                    |
| `src/rade_acq.h/c`                            | Acquisition and pilot detection                     |
| `src/rade_rx.h/c`                             | Receiver with sync state machine                    |
| `src/rade_api_nopy.c`                         | Python-free API implementation                      |
| `src/radae_tx_nopy.c` / `src/radae_rx_nopy.c` | Standalone executables                              |

The implementation uses built-in neural network weights (compiled from `rade_enc_data.c` and `rade_dec_data.c`), eliminating any need for Python or external model files at runtime.
