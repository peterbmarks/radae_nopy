# RADE V1 API Usage

This document describes the public C API of the Radio Autoencoder (RADE)
library as declared in [src/rade_api.h](src/rade_api.h). The API exposes a
single context that contains one transmitter (Tx) and one receiver (Rx).

All functions are declared `extern "C"` so the library can be used from both C
and C++. On Windows the symbols are exported/imported via `__declspec` and use
the `__stdcall` calling convention (handled automatically by the
`RADE_EXPORT` macro).

---

## Concepts and conventions

- **Complex samples** — IQ samples are passed as arrays of `RADE_COMP`, a
  single-precision complex type:

  ```c
  typedef struct {
    float real;
    float imag;
  } RADE_COMP;
  ```

- **Sample rates** — two fixed rates are defined:

  | Macro | Value | Meaning |
  | --- | --- | --- |
  | `RADE_MODEM_SAMPLE_RATE` | 8000 | modem (IQ) waveform sample rate |
  | `RADE_SPEECH_SAMPLE_RATE` | 16000 | speech sample rate |

- **Features** — RADE does not encapsulate the vocoder in this version. The
  caller supplies/consumes vocoder *feature* vectors (floats). `rade_tx()`
  accepts features and produces modem IQ samples; `rade_rx()` accepts modem IQ
  samples and produces features.

- **Soft-decision bits** — EOO (End-of-Over) bits are represented as floats in
  `+1.0` / `-1.0` form (i.e. *not* `1`/`0`).

- **Single context** — this version supports a single `struct rade` context,
  with exactly one Tx and one Rx inside it.

---

## Initialization and lifecycle

### `void rade_initialize(void)`

Must be called **before any other RADE function**. Initializes internal
library state.

### `void rade_finalize(void)`

Call when finished with RADE to release library-level resources.

### `struct rade *rade_open(char model_file[], int flags)`

Creates and returns a RADE context. `model_file` names the model weights file
to load. `flags` is a bitwise OR of:

| Flag | Value | Meaning |
| --- | --- | --- |
| `RADE_USE_C_ENCODER` | `0x1` | use the C encoder implementation |
| `RADE_USE_C_DECODER` | `0x2` | use the C decoder implementation |
| `RADE_FOFF_TEST` | `0x4` | frequency-offset test mode (developers only) |
| `RADE_VERBOSE_0` | `0x8` | reduce verbosity to "quiet" |

Returns a pointer to a `struct rade` context, used in all subsequent calls.

### `void rade_close(struct rade *r)`

Closes and frees a context previously returned by `rade_open()`.

### `int rade_version(void)`

Returns the API version, allowing callers to detect API changes.

### The `struct rade` context

```c
struct rade {
    int flags;
    int auxdata;
    int bottleneck;

    rade_tx_state tx;   /* Transmitter state */
    rade_rx_state rx;   /* Receiver state */
};
```

Treat this as an opaque handle through the API functions; the internal Tx/Rx
state structs are defined in [src/rade_tx.h](src/rade_tx.h) and
[src/rade_rx.h](src/rade_rx.h).

---

## Sizing helpers

Use these to size the arrays you pass to the Tx/Rx functions. They depend on
the context configuration, so call them after `rade_open()`.

| Function | Returns |
| --- | --- |
| `int rade_n_tx_out(struct rade *r)` | number of `RADE_COMP` samples produced per `rade_tx()` call |
| `int rade_n_tx_eoo_out(struct rade *r)` | number of `RADE_COMP` samples produced by `rade_tx_eoo()` |
| `int rade_nin_max(struct rade *r)` | maximum number of input samples ever required by `rade_rx()` (use for buffer allocation) |
| `int rade_n_features_in_out(struct rade *r)` | number of feature floats per Tx input frame / Rx output frame |
| `int rade_n_eoo_bits(struct rade *r)` | number of EOO soft-decision bits |

---

## Transmit (Tx)

### `int rade_tx(struct rade *r, RADE_COMP tx_out[], float features_in[])`

Encodes one frame of input `features_in[]` into modem IQ samples written to
`tx_out[]`.

- `features_in[]` must hold `rade_n_features_in_out(r)` floats.
- `tx_out[]` must have room for `rade_n_tx_out(r)` `RADE_COMP` samples.
- Returns the number of `RADE_COMP` samples written to `tx_out[]`.

> Note: the vocoder is not encapsulated by the API in this version — the caller
> is responsible for producing the feature vectors.

### `void rade_tx_set_eoo_bits(struct rade *r, float eoo_bits[])`

Sets the EOO bits that will be embedded in the End-of-Over frame. `eoo_bits[]`
must contain `rade_n_eoo_bits(r)` floats in `+1.0` / `-1.0` form (not `1`/`0`).

### Callsign helpers

```c
#define RADE_EOO_CALLSIGN_MAX 8
```

The EOO frame can carry a callsign. The maximum callsign length is
`RADE_EOO_CALLSIGN_MAX` (8) characters, not counting the null terminator
(8 chars × 7 bits = 56 bits, well within the 180 available EOO bits).

#### `void rade_tx_set_eoo_callsign(struct rade *r, const char *callsign)`

Encodes a callsign string into the EOO bits ready for transmission.

- `callsign` must be a null-terminated ASCII string of at most
  `RADE_EOO_CALLSIGN_MAX` characters.
- Only the first `RADE_EOO_CALLSIGN_MAX * 7` bits of the stored EOO array are
  overwritten; the remaining EOO bits are left unchanged.

#### `int rade_rx_get_eoo_callsign(const float *eoo_bits, int n_eoo_bits, char *callsign_out)`

Decodes a callsign from received EOO soft-decision bits.

- `eoo_bits`: array of `n_eoo_bits` floats in `+1` / `-1` form (as returned by
  `rade_rx()`).
- `callsign_out`: caller-supplied buffer of at least
  `RADE_EOO_CALLSIGN_MAX + 1` bytes.
- Returns the number of characters written, not counting the null terminator.

### `int rade_tx_eoo(struct rade *r, RADE_COMP tx_eoo_out[])`

Generates the final End-of-Over frame at the end of an over. `tx_eoo_out[]`
must have room for `rade_n_tx_eoo_out(r)` `RADE_COMP` samples. Returns the
number of `RADE_COMP` samples written.

---

## Receive (Rx)

### `int rade_nin(struct rade *r)`

Returns the number of input IQ samples to place into `rx_in[]` for the **next**
call to `rade_rx()`. This value can vary between calls (timing tracking), so
call `rade_nin()` before every `rade_rx()` call. Size your `rx_in[]` buffer
using `rade_nin_max()`.

### `int rade_rx(struct rade *r, float features_out[], int *has_eoo_out, float eoo_out[], RADE_COMP rx_in[])`

Processes `rade_nin(r)` input IQ samples from `rx_in[]`.

- `features_out[]`: receives decoded feature floats; size it with
  `rade_n_features_in_out(r)`.
- `has_eoo_out`: set by the function — non-zero means `eoo_out[]` contains
  valid End-of-Over soft-decision bits.
- `eoo_out[]`: receives EOO soft-decision bits from the QPSK symbols in
  `...IQIQI...` order; size it with `rade_n_eoo_bits(r)`.
- `rx_in[]`: input IQ samples (`rade_nin(r)` of them).

Returns non-zero if `features_out[]` contains valid output; the returned value
is the number of feature samples written to `features_out[]`.

### `int rade_sync(struct rade *r)`

Returns non-zero when the receiver is currently in sync.

### `float rade_freq_offset(struct rade *r)`

Returns the current frequency offset of the received signal. Valid only when
`rade_sync(r) != 0`.

### `int rade_snrdB_3k_est(struct rade *r)`

Returns the current SNR estimate in dB (3 kHz noise bandwidth) of the received
signal. Valid only when `rade_sync(r) != 0`.

### `void rade_set_disable_unsync(struct rade *r, float seconds)`

Test mode: disables the unsync mechanism for the given number of seconds
(`0` = disabled, i.e. normal unsync behavior).

---

## Typical usage

### Transmit loop

```c
rade_initialize();
struct rade *r = rade_open("model_file", RADE_USE_C_ENCODER | RADE_USE_C_DECODER);

int n_features = rade_n_features_in_out(r);
int n_tx       = rade_n_tx_out(r);

float     features_in[n_features];
RADE_COMP tx_out[n_tx];

/* optional: embed a callsign in the End-of-Over frame */
rade_tx_set_eoo_callsign(r, "N0CALL");

while (have_more_speech()) {
    fill_features(features_in);                 /* caller's vocoder */
    int nout = rade_tx(r, tx_out, features_in);
    send_iq(tx_out, nout);                       /* to the radio */
}

/* end of over */
RADE_COMP eoo_out[rade_n_tx_eoo_out(r)];
int neoo = rade_tx_eoo(r, eoo_out);
send_iq(eoo_out, neoo);

rade_close(r);
rade_finalize();
```

### Receive loop

```c
rade_initialize();
struct rade *r = rade_open("model_file", RADE_USE_C_ENCODER | RADE_USE_C_DECODER);

int n_features = rade_n_features_in_out(r);
int n_eoo      = rade_n_eoo_bits(r);

RADE_COMP rx_in[rade_nin_max(r)];
float     features_out[n_features];
float     eoo_out[n_eoo];

for (;;) {
    int nin = rade_nin(r);                       /* call before every rade_rx() */
    read_iq(rx_in, nin);                         /* from the radio */

    int has_eoo = 0;
    int nout = rade_rx(r, features_out, &has_eoo, eoo_out, rx_in);

    if (nout) {
        use_features(features_out, nout);        /* caller's vocoder */
    }
    if (has_eoo) {
        char callsign[RADE_EOO_CALLSIGN_MAX + 1];
        rade_rx_get_eoo_callsign(eoo_out, n_eoo, callsign);
        /* end of over detected; callsign decoded */
    }

    if (rade_sync(r)) {
        float foff = rade_freq_offset(r);
        int   snr  = rade_snrdB_3k_est(r);
        (void)foff; (void)snr;
    }
}

rade_close(r);
rade_finalize();
```

---

## Function quick reference

| Function | Purpose |
| --- | --- |
| `rade_initialize` | initialize library (call first) |
| `rade_finalize` | tear down library (call last) |
| `rade_open` | create a Tx/Rx context |
| `rade_close` | destroy a context |
| `rade_version` | query API version |
| `rade_n_tx_out` | Tx output samples per frame |
| `rade_n_tx_eoo_out` | Tx EOO output samples |
| `rade_nin_max` | max Rx input samples (buffer sizing) |
| `rade_n_features_in_out` | feature floats per frame |
| `rade_n_eoo_bits` | number of EOO bits |
| `rade_tx` | encode features → IQ samples |
| `rade_tx_set_eoo_bits` | set raw EOO bits |
| `rade_tx_set_eoo_callsign` | encode callsign into EOO bits |
| `rade_rx_get_eoo_callsign` | decode callsign from EOO bits |
| `rade_tx_eoo` | generate End-of-Over frame |
| `rade_nin` | input samples needed for next `rade_rx` |
| `rade_rx` | decode IQ samples → features (+ EOO) |
| `rade_sync` | is the Rx in sync? |
| `rade_freq_offset` | current Rx frequency offset |
| `rade_snrdB_3k_est` | current Rx SNR estimate (dB) |
| `rade_set_disable_unsync` | test mode: disable unsync |
