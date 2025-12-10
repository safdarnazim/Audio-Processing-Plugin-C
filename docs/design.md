# EP Audio Processing Plugin - Design Document

## Overview

The EP Audio Processing Plugin is a comprehensive C-based audio processing library designed to perform various digital signal processing (DSP) operations on WAV audio files. The system provides a menu-driven interface for applying professional audio effects including mixing, filtering, pitch shifting, denoising, and playback speed control.

## Architecture

### System Components

```
┌─────────────────────────────────────────────────────┐
│                  Main Menu System                    │
│              (User Interface Layer)                  │
└──────────────────┬──────────────────────────────────┘
                   │
      ┌────────────┴────────────┐
      │                         │
┌─────▼──────┐          ┌──────▼──────┐
│  WAV I/O   │          │   Plugin    │
│  Handler   │◄─────────┤   Manager   │
└────────────┘          └──────┬──────┘
                               │
        ┌──────────────────────┼──────────────────────┐
        │                      │                      │
   ┌────▼────┐          ┌──────▼──────┐       ┌──────▼──────┐
   │ Mixing  │          │  Filtering  │       │   Pitch     │
   │ Plugin  │          │   Plugin    │       │   Shift     │
   └─────────┘          └──────┬──────┘       └─────────────┘
                               │
                    ┌──────────┴──────────┐
                    │                     │
             ┌──────▼──────┐       ┌──────▼──────┐
             │  Low-Pass   │       │ High-Pass   │
             │   Filter    │       │   Filter    │
             └─────────────┘       └─────────────┘
```

## Core Data Structures

### WAVHeader Structure
```c
typedef struct {
    char riff_header[4];      // "RIFF"
    uint32_t wav_size;        // Total file size - 8 bytes
    char wave_header[4];      // "WAVE"
    char fmt_header[4];       // "fmt "
    uint32_t fmt_size;        // Format chunk size (16 for PCM)
    uint16_t audio_format;    // Audio format (1 = PCM)
    uint16_t no_of_channels;  // Number of channels (1=Mono, 2=Stereo)
    uint32_t sample_rate;     // Sample rate (Hz)
    uint32_t byte_rate;       // Byte rate
    uint16_t block_align;     // Block alignment
    uint16_t bit_depth;       // Bits per sample
    char data_header[4];      // "data"
    uint32_t data_bytes;      // Size of audio data
} WAVHeader;
```

## Plugin Modules

### 1. Audio Mixing Plugin (`mix_ep_plugin`)

**Purpose**: Combines two audio sources with adjustable gain control.

**Algorithm**:
- Loads two WAV files into memory
- Applies gain factors to each sample
- Mixes samples using the formula: `output = (gain1 × sample1) + (gain2 × sample2)`
- Implements clipping protection to prevent overflow

**Parameters**:
- `gain_sample1`: Gain multiplier for first audio source (default: 1.0)
- `gain_sample2`: Gain multiplier for second audio source (default: 1.5)

**Technical Details**:
- Uses 16-bit signed integer samples
- Implements saturation arithmetic (clips to ±32767)
- Processes samples in a single pass

### 2. Filtering Plugin (`filter_ep_plugin`)

#### Low-Pass Filter (LPF)

**Purpose**: Attenuates frequencies above a specified cutoff frequency.

**Algorithm**: FIR (Finite Impulse Response) filter using windowed-sinc method
- **Filter Design**: Sinc function with Hamming window
- **Cutoff Frequency**: 200 Hz (configurable)
- **Filter Size**: 101 taps

**Mathematical Foundation**:
```
h[n] = normalized_cutoff × sinc(π × (n - middle) × normalized_cutoff)
window[n] = 0.54 - 0.46 × cos(2π × n / (size - 1))
filter[n] = h[n] × window[n]
```

**Use Cases**:
- Removing high-frequency noise
- Bass extraction
- Rumble reduction

#### High-Pass Filter (HPF)

**Purpose**: Attenuates frequencies below a specified cutoff frequency.

**Algorithm**: First-order IIR (Infinite Impulse Response) filter
- **Cutoff Frequency**: 5000 Hz (configurable)
- **Filter Type**: RC High-Pass Filter

**Mathematical Foundation**:
```
RC = 1 / (2π × cutoff_frequency)
α = RC / (RC + dt)
output[n] = α × (output[n-1] + input[n] - input[n-1])
```

**Use Cases**:
- Removing low-frequency hum
- Isolating high-frequency content
- DC offset removal

### 3. Playback Speed Control (`plsp_ep_plugin`)

**Purpose**: Changes audio playback speed without affecting pitch.

**Algorithm**: Linear interpolation resampling
- Calculates new sample positions based on speed ratio
- Interpolates between adjacent samples for smooth transitions

**Speed Ratio**:
- `0.75` = 75% speed (slower)
- `1.0` = Normal speed
- `1.5` = 150% speed (faster)

**Mathematical Foundation**:
```
input_index = output_index / SPEED_RATIO
output[i] = (1 - fraction) × input[index1] + fraction × input[index2]
```

**Note**: This implementation changes duration without pitch correction. For pitch-preserving speed changes, additional phase vocoder techniques would be needed.

### 4. Pitch Shifting Plugin (`pitchshift_ep_plugin`)

**Purpose**: Changes the pitch of audio by a specified number of semitones.

**Algorithm**: Resampling with frequency scaling
- Uses equal temperament tuning system
- Applies linear interpolation for sample rate conversion

**Mathematical Foundation**:
```
pitch_factor = 2^(semitones / 12)
new_length = original_length / pitch_factor
```

**Default Configuration**:
- Shift: +8 semitones (up a major sixth)

**Semitone Examples**:
- `-12`: Down one octave
- `-7`: Down a perfect fifth
- `0`: No change
- `+7`: Up a perfect fifth
- `+12`: Up one octave

**Limitation**: This simple implementation changes both pitch and duration. Professional pitch shifting requires time-stretching algorithms (e.g., PSOLA, Phase Vocoder).

### 5. Denoising Plugin (`denoise_ep_plugin`)

**Purpose**: Reduces unwanted noise in audio recordings.

**Algorithm**: Low-pass filtering approach
- Applies FIR filter with 3000 Hz cutoff
- Removes high-frequency noise components
- Preserves speech and music content in lower frequencies

**Technical Details**:
- Filter Size: 101 taps
- Cutoff Frequency: 3000 Hz
- Method: Windowed-sinc with Hamming window

**Best For**:
- Hiss reduction
- White noise removal
- High-frequency interference

**Limitations**: Simple spectral filtering; advanced denoising would use spectral subtraction or Wiener filtering.

### 6. Audio Information Reader (`get_info`)

**Purpose**: Displays technical specifications of WAV files.

**Information Provided**:
- File size (bytes)
- Sample rate (Hz)
- Audio format (PCM = 1)
- Number of channels
- Bit depth
- Audio data size (bytes)

## Configuration Parameters

### Global Settings
```c
#define SAMPLE_RATE 48000           // Default sample rate: 48 kHz
#define PI 3.14159265358979323846   // Mathematical constant
```

### Input/Output Paths
```c
#define input_sample1 "audio_files/bgm.wav"
#define input_sample2 "audio_files/vocals.wav"
#define filter_in "audio_files/vocals.wav"
#define noise_input "audio_files/noise_audio.wav"
#define pitch_input "audio_files/vocals.wav"
```

### Processing Parameters
```c
#define gain_sample1 1.0          // Mixing gain for input 1
#define gain_sample2 1.5          // Mixing gain for input 2
#define SPEED_RATIO 0.75          // Playback speed multiplier
#define SEMITONE_FACTOR(s) (pow(2.0, (s) / 12.0))  // Pitch shift calculation
```

## Signal Processing Techniques

### 1. Linear Interpolation
Used in pitch shifting and playback speed control for smooth sample transitions.

**Advantages**:
- Simple implementation
- Low computational cost
- Smooth transitions

**Disadvantages**:
- May introduce high-frequency artifacts
- Not ideal for large pitch shifts

### 2. Windowed-Sinc Filter Design
Used in low-pass filtering and denoising.

**Advantages**:
- Excellent frequency response
- Sharp cutoff characteristics
- Adjustable filter order

**Disadvantages**:
- Requires convolution (computationally intensive)
- Introduces latency

### 3. IIR Filtering
Used in high-pass filtering.

**Advantages**:
- Efficient (requires less memory)
- Real-time processing friendly
- Simple implementation

**Disadvantages**:
- Can be unstable if not designed carefully
- Phase distortion

## Memory Management

### Buffer Allocation Strategy
- Dynamic allocation using `malloc()`
- Proper cleanup with `free()`
- Error handling for allocation failures

### Sample Processing
- Processes entire audio files in memory
- Suitable for files up to available RAM
- For large files, consider streaming architecture

## Error Handling

The system implements comprehensive error checking:
- File I/O validation
- Memory allocation verification
- Boundary condition checks
- User input validation

## File Format Support

### Supported Format
- **WAV (WAVE)**: Uncompressed PCM audio
  - 16-bit signed integer samples
  - Mono or Stereo
  - 48 kHz sample rate (default)

### Unsupported Formats
- MP3, AAC, FLAC, OGG (compressed formats)
- Non-PCM WAV variants
- 24-bit or 32-bit depth audio

## Performance Considerations

### Computational Complexity
- **Mixing**: O(n) - Linear time
- **Filtering**: O(n × m) - Linear in samples, linear in filter size
- **Pitch Shift**: O(n) - Linear with interpolation
- **Playback Speed**: O(n) - Linear with interpolation

### Optimization Opportunities
1. **SIMD Instructions**: Vectorize sample processing
2. **Multi-threading**: Process channels in parallel
3. **Block Processing**: Process in chunks for better cache utilization
4. **Fixed-Point Math**: Replace floating-point in filters

## Usage Workflow

```
1. Program Start
   ↓
2. Display Menu
   ↓
3. User Selects Operation
   ↓
4. Load Audio File(s)
   ↓
5. Apply DSP Algorithm
   ↓
6. Write Output File
   ↓
7. Return to Menu (Loop)
```

## Future Enhancement Possibilities

### Planned Features
1. **Real-time Processing**: Stream-based processing for large files
2. **Advanced Pitch Shifting**: Phase vocoder implementation
3. **Equalizer**: Multi-band parametric EQ
4. **Compressor/Limiter**: Dynamic range control
5. **Reverb**: Convolution-based room simulation
6. **Batch Processing**: Process multiple files automatically

### Additional File Format Support
- MP3 encoding/decoding (via external library)
- FLAC support for lossless compression
- Multi-file format conversion

### GUI Development
- Cross-platform desktop application
- Real-time waveform visualization
- Drag-and-drop file support

## Dependencies

### Standard Libraries
- `stdio.h` - File I/O operations
- `stdlib.h` - Memory management
- `stdint.h` - Fixed-width integer types
- `math.h` - Mathematical functions (sin, cos, pow)

### External Dependencies
None - Pure C implementation with standard libraries only.

## Build Instructions

### Compilation
```bash
gcc -o audio_plugin DSP.c -lm
```

The `-lm` flag links the math library for mathematical functions.

### Running
```bash
./audio_plugin
```

## Testing Recommendations

### Unit Testing
1. Test each plugin independently
2. Verify WAV header parsing
3. Validate sample-level operations
4. Check boundary conditions

### Integration Testing
1. Chain multiple effects
2. Test with various audio formats (mono/stereo, different sample rates)
3. Validate file I/O with edge cases
4. Memory leak detection (valgrind)

### Audio Quality Testing
1. Compare input/output spectrograms
2. Measure SNR (Signal-to-Noise Ratio)
3. Listen tests for artifacts
4. Frequency response analysis

## License & Credits

This audio processing plugin demonstrates fundamental DSP concepts in C. It's designed for educational purposes and can be extended for production use with additional error handling and optimization.

---

**Version**: 1.0  
**Last Updated**: December 2025  
**Language**: C (C99 standard)  
**Platform**: Cross-platform (Linux, Windows, macOS)