# Audio Processing Plugin C

This project implements a collection of audio processing features in C for WAV files. It provides mixing, playback speed control, pitch shifting, denoising, low-pass filtering, high-pass filtering, and WAV file information extraction. All operations are implemented using custom DSP logic without third-party libraries.

---

## Features

1. **Mixing** of two WAV audio files with independent gain control
2. **Playback speed modification** using interpolation-based resampling
3. **Pitch shifting** in semitone steps using frequency scaling
4. **Denoising** using a FIR low-pass filter
5. **Low-pass filtering** using a windowed-sinc FIR filter
6. **High-pass filtering** using a simple IIR filter
7. **WAV metadata display** (sample rate, channels, bit depth, etc.)

---

## Project Structure

```
Audio-Processing-Plugin-C/
│
├── src/
│   └── audio_plugin.c
│
├── audio_files/
│   ├── bgm.wav
│   ├── vocals.wav
│   ├── output_denoised2.wav
│   └── README.txt
│
├── docs/
│   ├── design.md
│   └── screenshots/
│
├── Makefile
├── .gitignore
├── LICENSE
└── README.md
```

---

## Building and Running

### Build

```bash
make
```

### Run

```bash
make run
```

### Manual Compilation

```bash
gcc src/audio_plugin.c -o audio_plugin -lm
./audio_plugin
```

---

## Audio Files

The `audio_files/` folder contains sample WAV files for testing.

Generated output files include:

- `output.wav`
- `output_lpf.wav`
- `output_hpf.wav`
- `output_denoise.wav`
- `output_pitch_shift.wav`
- `output_psb.wav`

You may replace these with your own WAV files (must be PCM format).

---

## Documentation

See `docs/design.md` for internal DSP explanation, filter design, and processing flow.

---

## Usage

The program provides an interactive menu with the following options:

1. **Mix two audio files** - Combine two WAV files with adjustable gain
2. **Change playback speed** - Speed up or slow down audio
3. **Pitch shift** - Change pitch by semitones
4. **Denoise audio** - Remove noise using low-pass filtering
5. **Apply low-pass filter** - Filter out high frequencies
6. **Apply high-pass filter** - Filter out low frequencies
7. **Display WAV info** - Show metadata of a WAV file
8. **Exit** - Quit the program

---

## Requirements

- GCC compiler
- Standard C library
- Math library (`-lm` flag)
- PCM format WAV files

---

## Technical Details

### Audio Processing Algorithms

- **Mixing**: Linear mixing with gain control
- **Speed Control**: Linear interpolation resampling
- **Pitch Shifting**: Frequency domain scaling
- **Denoising**: FIR low-pass filter
- **Low-Pass Filter**: Windowed-sinc FIR filter
- **High-Pass Filter**: Simple IIR filter

All DSP operations are implemented from scratch without external audio processing libraries.

---

## License

This project is licensed under the MIT License.

---

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

---

## Author

Created as an educational project demonstrating audio DSP concepts in C.

---

## Notes

- All audio files must be in WAV format with PCM encoding
- Ensure sufficient disk space for output files
- Processing time depends on file size and selected operation
- Backup original files before processing

---

## Future Enhancements

Potential features for future versions:

- Multi-band EQ
- Compression/limiting
- Reverb effects
- Real-time processing
- GUI interface
- Support for additional audio formats

---

## Troubleshooting

**Issue**: Compilation errors
- **Solution**: Ensure GCC is installed and `-lm` flag is included

**Issue**: "File not found" errors
- **Solution**: Verify audio files exist in the `audio_files/` directory

**Issue**: Distorted output
- **Solution**: Check gain levels and input file format

**Issue**: Unexpected behavior
- **Solution**: Ensure input files are valid PCM WAV files

---

## References

For more information on digital signal processing and audio programming, refer to `docs/design.md`.
