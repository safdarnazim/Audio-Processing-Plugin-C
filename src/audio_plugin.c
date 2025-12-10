#include<stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

//MIXING
#define input_sample1 "audio_files/bgm.wav" //INPUT AUDIO 1
#define gain_sample1 1.0 //GAIN OF INPUT 1
#define input_sample2 "audio_files/vocals.wav" //INPUT AUDIO 2
#define gain_sample2 1.5 //GAIN OF INPUT 2

//PLAYBACK CONTROL
#define SPEED_RATIO 0.75 //ADJUST THIS RATIO  TO CHANGE PLAYBACK SPEED

//DENOISE
#define SAMPLE_RATE 48000
#define PI 3.14159265358979323846
#define noise_input "audio_files/onoise_audio.wav"

//FILTER INPUT
#define filter_in "audio_files/vocals.wav"

//PITCH SHIFTER
#define SEMITONE_FACTOR(semitones) (pow(2.0, (semitones) / 12.0))
#define pitch_input "audio_files/vocals.wav"

//#define action 1//

void mix_ep_plugin();
void filter_ep_plugin();
void pitchshift_ep_plugin();
void plsp_ep_plugin();
void denoise_ep_plugin();
void get_info();

int main()
{
    int action;
    int flag=1;

    printf("\t\t\t\t\tEP AUDIO PROCESsING PLUGIN\n\n");
    printf("Choose the action you need to apply on you WAV audio file/s\n\n");
    printf("1.Mixing\n");
    printf("2.Playback speed control\n");
    printf("3.Pitch Shift\n");
    printf("4.Denoise\n");
    printf("5.Filtering\n");
    printf("6.Get input Audio File Info\n");

    while(flag==1)
    {
        printf("\n\n\n");
        printf("Enter the code number for action: ");
        scanf("%d",&action);

        switch(action)
        {
            case 1:mix_ep_plugin(input_sample1,input_sample2,"audio_files/output.wav");
                   printf("Succesfully mixed!! ENJOY\n\n");
                   break;

            case 2:plsp_ep_plugin(input_sample1);
            break;

            case 3:pitchshift_ep_plugin(pitch_input);
            break;

            case 4:denoise_ep_plugin(noise_input);
            break;

            case 5:filter_ep_plugin(filter_in);
            break;

            case 6:get_info(input_sample2);
            break;

            default:printf("ERROR CODE\n\n");
            break;
        }

        printf("______________x________________x_________________x_________________");
    }

}

//MIX EP PLUGIN start
    typedef struct{
        char riff_header[4];
        uint32_t wav_size; //Total file size-8
        char wave_header[4];
        char fmt_header[4];
        uint32_t fmt_size;
        uint16_t audio_format;
        uint16_t no_of_channels;
        uint32_t sample_rate;
        uint32_t byte_rate;
        uint16_t block_align;
        uint16_t bit_depth;
        char data_header[4];
        uint32_t data_bytes;
    }WAVHeader;

    void wav_header_read(FILE *file,WAVHeader *header)
    {
        fread(header,sizeof(WAVHeader),1,file);
    }

    void wav_header_write(FILE *file,WAVHeader *header)
    {
        fwrite(header,sizeof(WAVHeader),1,file);
    }

   int16_t mix_samples(int16_t sample1, int16_t sample2) {
    int32_t mixed_sample = (gain_sample1 * sample1) + (gain_sample2 * sample2);
    if (mixed_sample > 32767) mixed_sample = 32767;
    if (mixed_sample < -32768) mixed_sample = -32768;
    return (int16_t)mixed_sample;
}

void mix_ep_plugin(const char* input_file1, const char* input_file2, const char* output_file_path) {
    // Open input files
    FILE *file1 = fopen(input_file1, "rb");
    FILE *file2 = fopen(input_file2, "rb");

    if (!file1 || !file2) {
        printf("Error opening input files.\n");
        if (file1) fclose(file1);
        if (file2) fclose(file2);
        return;
    }

    // Read headers
    WAVHeader header1, header2;
    wav_header_read(file1, &header1);
    wav_header_read(file2, &header2);

    // Calculate number of samples (use the smaller of the two)
    int num_samples1 = header1.data_bytes / (header1.bit_depth / 8);
    int num_samples2 = header2.data_bytes / (header2.bit_depth / 8);
    int num_samples = (num_samples1 < num_samples2) ? num_samples1 : num_samples2;

    // Allocate memory for samples
    int16_t *input1_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));
    int16_t *input2_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));
    int16_t *output_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));

    if (!input1_samples || !input2_samples || !output_samples) {
        printf("Memory allocation failed.\n");
        fclose(file1);
        fclose(file2);
        free(input1_samples);
        free(input2_samples);
        free(output_samples);
        return;
    }

    // Read audio data
    fread(input1_samples, sizeof(int16_t), num_samples, file1);
    fread(input2_samples, sizeof(int16_t), num_samples, file2);
    fclose(file1);
    fclose(file2);

    // Mix samples
    for (int i = 0; i < num_samples; i++) {
        output_samples[i] = mix_samples(input1_samples[i], input2_samples[i]);
    }

    // Write output file
    FILE *output_file = fopen(output_file_path, "wb");
    if (!output_file) {
        printf("Error opening output file.\n");
        free(input1_samples);
        free(input2_samples);
        free(output_samples);
        return;
    }

    // Update header for output
    header1.data_bytes = num_samples * sizeof(int16_t);
    header1.wav_size = header1.data_bytes + sizeof(WAVHeader) - 8;

    wav_header_write(output_file, &header1);
    fwrite(output_samples, sizeof(int16_t), num_samples, output_file);

    fclose(output_file);
    free(input1_samples);
    free(input2_samples);
    free(output_samples);
}
//MIX EP PLUGIN ENDS

//GET INFO START
void get_info(input_file)
{
    printf("\tSAMPLE INFORMATIONS \n");
    FILE *file1 = fopen(input_file, "rb");
    WAVHeader header1;
    wav_header_read(file1, &header1);
    printf("Overall Size=%d bytes \nSample Rate= %d Hz \nAudio Format =%d \nNo of Channels =%d \nBit Depth =%d \nData Size=%d bytes\n\n",
           header1.wav_size,
           header1.sample_rate,
           header1.audio_format,
           header1.no_of_channels,
           header1.bit_depth,
           header1.data_bytes);
}
//GET INFO END

//Filter STARTs
void filter_ep_plugin(input_file)
{
    int filter_action;
    printf("\n1.Low Pass Filtering\n2.High Pass Filtering\nEnter your Action Number:");
    scanf("%d",&filter_action);
    if(filter_action==1)
    {
        lpfilter_ep_plugin(input_file);
    }

    else if(filter_action==2)
    {
        hpfilter_ep_plugin(input_file);
    }
    else
        printf("\nError CODE\n");
}

void design_lpf(float* filter, int filter_size, float cutoff_freq, int sample_rate) {
    float norm_cutoff = 2.0 * cutoff_freq / sample_rate;  // Normalized cutoff frequency
    int middle = filter_size / 2;

    for (int i = 0; i < filter_size; i++) {
        if (i == middle) {
            filter[i] = norm_cutoff;
        } else {
            float sinc_val = sin(PI * (i - middle) * norm_cutoff) / (PI * (i - middle));
            filter[i] = sinc_val;
        }
        // Apply a Hamming window to reduce side lobes
        filter[i] *= (0.54 - 0.46 * cos(2 * PI * i / (filter_size - 1)));
    }
}

void apply_lpf(int16_t* input_samples, int num_samples, int16_t* output_samples, float* filter, int filter_size) {
    for (int i = 0; i < num_samples; i++) {
        float sum = 0.0;
        for (int j = 0; j < filter_size; j++) {
            if (i - j >= 0) {
                sum += input_samples[i - j] * filter[j];
            }
        }
        output_samples[i] = (int16_t)sum;
    }
}

void lpfilter_ep_plugin(input_files)
{
    FILE *input_file = fopen(input_files, "rb");
    if (!input_file) {
        printf("Error opening input file.\n");
        return 1;
    }

    // Read WAV header (assuming WAVHeader structure exists)
    WAVHeader header;
    fread(&header, sizeof(WAVHeader), 1, input_file);

    // Read audio samples
    int num_samples = header.data_bytes / (header.bit_depth / 8);
    int16_t *input_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));
    fread(input_samples, sizeof(int16_t), num_samples, input_file);
    fclose(input_file);

    // Set filter parameters
    int filter_size = 101;  // Size of the filter (higher = more precision)
    float cutoff_freq = 200.0;  // Cutoff frequency (Hz)
    float *lpf_filter = (float*)malloc(filter_size * sizeof(float));

    // Design the low-pass filter
    design_lpf(lpf_filter, filter_size, cutoff_freq, 48000);  // 48 kHz sample rate

    // Apply the filter
    int16_t *output_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));
    apply_lpf(input_samples, num_samples, output_samples, lpf_filter, filter_size);

    // Write the filtered output to a new WAV file
    FILE *output_file = fopen("audio_files/output_lpf.wav", "wb");
    if (!output_file) {
        printf("Error opening output file.\n");
        free(input_samples);
        free(output_samples);
        free(lpf_filter);
        return 1;
    }

    // Update the WAV header for the new file
    fwrite(&header, sizeof(WAVHeader), 1, output_file);
    fwrite(output_samples, sizeof(int16_t), num_samples, output_file);

    // Clean up
    fclose(output_file);
    free(input_samples);
    free(output_samples);
    free(lpf_filter);

    printf("Low-pass filtering completed successfully.\n");
}

void high_pass_filter(int16_t* input, int16_t* output, int num_samples, float cutoff_freq, float sample_rate) {
    float dt = 1.0 / sample_rate;
    float RC = 1.0 / (2 * PI * cutoff_freq);
    float alpha = RC / (RC + dt);

    // Initialize previous input and output
    float prev_input = input[0];
    float prev_output = input[0];

    for (int i = 1; i < num_samples; i++) {
        // Apply high-pass filter equation with the alpha coefficient
        output[i] = (int16_t)(alpha * (prev_output + input[i] - prev_input));

        // Update previous input and output
        prev_output = output[i];
        prev_input = input[i];
    }
}

void hpfilter_ep_plugin(input_files)
{
    FILE *input_file = fopen(input_files, "rb");
    if (!input_file) {
        printf("Error opening input file.\n");
        return 1;
    }

    // Read WAV header
    WAVHeader header;
    fread(&header, sizeof(WAVHeader), 1, input_file);

    // Calculate number of input samples
    int num_samples = header.data_bytes / (header.bit_depth / 8);
    int16_t *input_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));
    int16_t *output_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));

    if (input_samples == NULL || output_samples == NULL) {
        printf("Memory allocation failed.\n");
        fclose(input_file);
        return 1;
    }

    // Read input samples
    fread(input_samples, sizeof(int16_t), num_samples, input_file);
    fclose(input_file);

    // Define cutoff frequency and sample rate
    float cutoff_frequency = 5000.0;  // Cutoff frequency for HPF in Hz
    float sample_rate = 48000.0;      // Assuming 48kHz sample rate

    // Apply the high-pass filter
    high_pass_filter(input_samples, output_samples, num_samples, cutoff_frequency, sample_rate);

    // Write the filtered output to a new WAV file
    FILE *output_file = fopen("audio_files/output_hpf.wav", "wb");
    if (!output_file) {
        printf("Error opening output file.\n");
        free(input_samples);
        free(output_samples);
        return 1;
    }

    // Update the WAV header for the new file
    fwrite(&header, sizeof(WAVHeader), 1, output_file);

    // Write filtered samples to the output file
    fwrite(output_samples, sizeof(int16_t), num_samples, output_file);

    fclose(output_file);

    // Free memory
    free(input_samples);
    free(output_samples);

    printf("High-pass filter applied successfully.\n");
    return 0;
}
//FILTER ENDS


//PLAYBACK SPEED START
void plsp_ep_plugin()
{
    void linear_interpolation(int16_t* input_samples, int num_input_samples, int16_t* output_samples, int num_output_samples) {
    for (int i = 0; i < num_output_samples; i++) {
        float input_index = i / SPEED_RATIO;
        int index1 = (int)input_index;  // Floor of input_index
        int index2 = index1 + 1;        // Next sample

        // Boundary check for index2
        if (index2 >= num_input_samples) {
            output_samples[i] = input_samples[index1];
        } else {
            float fraction = input_index - index1;
            output_samples[i] = (int16_t)((1.0 - fraction) * input_samples[index1] + fraction * input_samples[index2]);
        }
    }
}
    FILE *input_file = fopen("audio_files/vocals.wav", "rb");
    if (!input_file) {
        printf("Error opening input file.\n");
        return 1;
    }

    // Read WAV header
    WAVHeader header;
    fread(&header, sizeof(WAVHeader), 1, input_file);

    // Calculate the number of input samples
    int num_input_samples = header.data_bytes / (header.bit_depth / 8);
    int16_t *input_samples = (int16_t*)malloc(num_input_samples * sizeof(int16_t));

    if (input_samples == NULL) {
        printf("Memory allocation failed for input samples.\n");
        fclose(input_file);
        return 1;
    }

    // Read the audio data
    fread(input_samples, sizeof(int16_t), num_input_samples, input_file);
    fclose(input_file);

    // Calculate number of output samples based on SPEED_RATIO
    int num_output_samples = (int)(num_input_samples * SPEED_RATIO);
    int16_t *output_samples = (int16_t*)malloc(num_output_samples * sizeof(int16_t));

    if (output_samples == NULL) {
        printf("Memory allocation failed for output samples.\n");
        free(input_samples);
        return 1;
    }

    // Perform linear interpolation to resample
    linear_interpolation(input_samples, num_input_samples, output_samples, num_output_samples);

    // Write output to a new WAV file
    FILE *output_file = fopen("audio_files/output_psb.wav", "wb");
    if (!output_file) {
        printf("Error opening output file.\n");
        free(input_samples);
        free(output_samples);
        return 1;
    }

    // Update the WAV header for the new file (adjust data_size and overall_size)
    header.data_bytes = num_output_samples * (header.bit_depth / 8);
    header.wav_size = header.data_bytes + sizeof(WAVHeader) - 8;

    // Write updated header to output file
    fwrite(&header, sizeof(WAVHeader), 1, output_file);

    // Write the resampled audio data
    fwrite(output_samples, sizeof(int16_t), num_output_samples, output_file);
    fclose(output_file);

    // Free memory
    free(input_samples);
    free(output_samples);

    printf("Playback speed adjusted SUCCESSFULLY'.\n");
    return 0;
}
//PLAYBACK SPEED ENDS

//DENOISE START
void design_denoise(float* filter, int filter_size, float cutoff_freq, int sample_rate) {
    float norm_cutoff = 2.0 * cutoff_freq / sample_rate;  // Normalized cutoff frequency
    int middle = filter_size / 2;

    for (int i = 0; i < filter_size; i++) {
        if (i == middle) {
            filter[i] = norm_cutoff;
        } else {
            float sinc_val = sin(PI * (i - middle) * norm_cutoff) / (PI * (i - middle));
            filter[i] = sinc_val;
        }
        // Apply a Hamming window to reduce side lobes
        filter[i] *= (0.54 - 0.46 * cos(2 * PI * i / (filter_size - 1)));
    }
}

void apply_denoise(int16_t* input_samples, int num_samples, int16_t* output_samples, float* filter, int filter_size) {
    for (int i = 0; i < num_samples; i++) {
        float sum = 0.0;
        for (int j = 0; j < filter_size; j++) {
            if (i - j >= 0) {
                sum += input_samples[i - j] * filter[j];
            }
        }
        output_samples[i] = (int16_t)sum;
    }
}

void denoise_ep_plugin(input_files)
{
    FILE *input_file = fopen(input_files, "rb");
    if (!input_file) {
        printf("Error opening input file.\n");
        return 1;
    }

    // Read WAV header (assuming WAVHeader structure exists)
    WAVHeader header;
    fread(&header, sizeof(WAVHeader), 1, input_file);

    // Read audio samples
    int num_samples = header.data_bytes / (header.bit_depth / 8);
    int16_t *input_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));
    fread(input_samples, sizeof(int16_t), num_samples, input_file);
    fclose(input_file);

    // Set filter parameters
    int filter_size = 101;  // Size of the filter (higher = more precision)
    float cutoff_freq = 3000.0;  // Cutoff frequency (Hz)
    float *lpf_filter = (float*)malloc(filter_size * sizeof(float));

    // Design the low-pass filter
    design_denoise(lpf_filter, filter_size, cutoff_freq, 48000);  // 48 kHz sample rate

    // Apply the filter
    int16_t *output_samples = (int16_t*)malloc(num_samples * sizeof(int16_t));
    apply_denoise(input_samples, num_samples, output_samples, lpf_filter, filter_size);

    // Write the filtered output to a new WAV file
    FILE *output_file = fopen("audio_files/output_denoise.wav", "wb");
    if (!output_file) {
        printf("Error opening output file.\n");
        free(input_samples);
        free(output_samples);
        free(lpf_filter);
        return 1;
    }

    // Update the WAV header for the new file
    fwrite(&header, sizeof(WAVHeader), 1, output_file);
    fwrite(output_samples, sizeof(int16_t), num_samples, output_file);

    // Clean up
    fclose(output_file);
    free(input_samples);
    free(output_samples);
    free(lpf_filter);

    printf("Denoise completed successfully.\n");
}
//DENOISE ENDS

//PTICH SHIFT START
void pitchshift_ep_plugin(input_files)
{
    void linear_interpolation(int16_t* input_samples, int num_input_samples, int16_t* output_samples, int num_output_samples, float pitch_factor) {
        for (int i = 0; i < num_output_samples; i++) {
            float input_index = i / pitch_factor;
            int index1 = (int)input_index;  // Floor of input_index
            int index2 = index1 + 1;        // Next sample

            // Boundary check for index2
            if (index2 >= num_input_samples) {
                output_samples[i] = input_samples[index1];
            } else {
                float fraction = input_index - index1;
                output_samples[i] = (int16_t)((1.0 - fraction) * input_samples[index1] + fraction * input_samples[index2]);
            }
        }
    }

    int pitch_shift(const char* input_file_path, const char* output_file_path, int semitone_shift) {
        // Open input file
        FILE* input_file = fopen(input_file_path, "rb");
        if (!input_file) {
            printf("Error opening input file.\n");
            return 1;
        }

        // Read WAV header
        WAVHeader header;
        fread(&header, sizeof(WAVHeader), 1, input_file);

        // Calculate the number of input samples
        int num_input_samples = header.data_bytes / (header.bit_depth / 8);
        int16_t* input_samples = (int16_t*)malloc(num_input_samples * sizeof(int16_t));

        if (input_samples == NULL) {
            printf("Memory allocation failed for input samples.\n");
            fclose(input_file);
            return 1;
        }

        // Read the audio data
        fread(input_samples, sizeof(int16_t), num_input_samples, input_file);
        fclose(input_file);

        // Calculate the pitch factor (scaling factor for the frequency)
        float pitch_factor = SEMITONE_FACTOR(semitone_shift);

        // Calculate number of output samples based on the pitch factor
        int num_output_samples = (int)(num_input_samples / pitch_factor);
        int16_t* output_samples = (int16_t*)malloc(num_output_samples * sizeof(int16_t));

        if (output_samples == NULL) {
            printf("Memory allocation failed for output samples.\n");
            free(input_samples);
            return 1;
        }

        // Perform linear interpolation for pitch shifting
        linear_interpolation(input_samples, num_input_samples, output_samples, num_output_samples, pitch_factor);

        // Open output file
        FILE* output_file = fopen(output_file_path, "wb");
        if (!output_file) {
            printf("Error opening output file.\n");
            free(input_samples);
            free(output_samples);
            return 1;
        }

        // Update the WAV header for the new file (adjust data_bytes and wav_size)
        header.data_bytes = num_output_samples * (header.bit_depth / 8);
        header.wav_size = header.data_bytes + sizeof(WAVHeader) - 8;

        // Write updated header to output file
        fwrite(&header, sizeof(WAVHeader), 1, output_file);

        // Write the pitch-shifted audio data
        fwrite(output_samples, sizeof(int16_t), num_output_samples, output_file);
        fclose(output_file);

        // Free memory
        free(input_samples);
        free(output_samples);

        printf("Pitch shifted successfully.\n");
        return 0;
    }

    int semitone_shift = 8; // Transpose up by 2 semitones
    pitch_shift(input_files, "audio_files/output_pitch_shift.wav", semitone_shift);
    return 0;
}

//PITCH SHIFT ENDs
