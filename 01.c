#include <math.h>
#include <string.h>

#define FS 360.0f          // Sampling frequency in Hz
#define WINDOW_SIZE 30     // Integration window size
#define BUFFER_SIZE 900    // Number of ECG samples
#define MAX_PEAKS 100      // Maximum number of peaks

// Buffers
float ecg_signal[BUFFER_SIZE];        // Input ECG signal
float y_lp[BUFFER_SIZE];              // Low-pass output
float y_hp[BUFFER_SIZE];              // High-pass output
float y_der[BUFFER_SIZE];             // Derivative output
float y_sq[BUFFER_SIZE];              // Squared output
float y_int[BUFFER_SIZE];             // Integrated signal
int r_peak_binary[BUFFER_SIZE];       // Binary array with R-peaks
int r_peaks[MAX_PEAKS];               // Indices of detected R-peaks
float heart_rate = 0.0f;
int num_peaks = 0;

// Filter and processing stages
void low_pass_filter(float *x, float *y, int N) {
    memset(y, 0, N * sizeof(float));
    for (int n = 12; n < N; n++) {
        y[n] = 2 * y[n - 1] - y[n - 2] + x[n] - 2 * x[n - 6] + x[n - 12];
    }
}

void high_pass_filter(float *x, float *y, int N) {
    memset(y, 0, N * sizeof(float));
    for (int n = 32; n < N; n++) {
        y[n] = y[n - 1] - (1.0f / 32.0f) * x[n] + x[n - 16] - x[n - 17] + (1.0f / 32.0f) * x[n - 32];
    }
}

void derivative(float *x, float *y, int N) {
    memset(y, 0, N * sizeof(float));
    for (int n = 4; n < N; n++) {
        y[n] = (1.0f / 8.0f) * (2 * x[n] + x[n - 1] - x[n - 3] - 2 * x[n - 4]);
    }
}

void square(float *x, float *y, int N) {
    for (int n = 0; n < N; n++) {
        y[n] = x[n] * x[n];
    }
}

void integration(float *x, float *y, int N) {
    memset(y, 0, N * sizeof(float));
    for (int n = WINDOW_SIZE - 1; n < N; n++) {
        float sum = 0.0f;
        for (int i = 0; i < WINDOW_SIZE; i++) {
            sum += x[n - i];
        }
        y[n] = sum / WINDOW_SIZE;
    }
}

void peak_detection(float *x, int N, float fs, int *r_peaks, int *r_peak_binary, int *num_peaks, float *heart_rate) {
    float SPKI = 0.0f, NPKI = 0.0f;
    float THRESHOLD1, THRESHOLD2;
    int last_peak = 0;
    *num_peaks = 0;

    memset(r_peak_binary, 0, N * sizeof(int));

    float max_val = x[0];
    for (int i = 1; i < N; i++) {
        if (x[i] > max_val) max_val = x[i];
    }
    SPKI = 0.25f * max_val;
    NPKI = 0.125f * max_val;
    THRESHOLD1 = NPKI + 0.25f * (SPKI - NPKI);
    THRESHOLD2 = 0.5f * THRESHOLD1;

    for (int n = 1; n < N - 1; n++) {
        if (x[n] > x[n - 1] && x[n] > x[n + 1]) {
            float PEAKI = x[n];
            if (PEAKI > THRESHOLD1) {
                if (*num_peaks == 0 || (n - last_peak) > (int)(0.6f * fs)) {
                    if (*num_peaks < MAX_PEAKS) {
                        r_peaks[*num_peaks] = n;
                        r_peak_binary[n] = 1;
                        (*num_peaks)++;
                        last_peak = n;
                        SPKI = 0.125f * PEAKI + 0.875f * SPKI;
                    }
                }
            } else {
                NPKI = 0.125f * PEAKI + 0.875f * NPKI;
            }
            THRESHOLD1 = NPKI + 0.25f * (SPKI - NPKI);
            THRESHOLD2 = 0.5f * THRESHOLD1;
        }
    }

    if (*num_peaks > 1) {
        float rr_sum = 0.0f;
        for (int i = 1; i < *num_peaks; i++) {
            float rr = (float)(r_peaks[i] - r_peaks[i - 1]) / fs;
            rr_sum += 60.0f / rr;
        }
        *heart_rate = rr_sum / (*num_peaks - 1);
    } else {
        *heart_rate = 0.0f;
    }
}

int main() {
    // Manually load ECG samples here (or load them from array if importing from text isn't supported)
    // Example for quick testing:
    //for (int i = 0; i < BUFFER_SIZE; i++) {
    //    ecg_signal[i] = sinf(2 * M_PI * 1.0f * i / FS);  // Replace with actual ECG data
    //}
    FILE *fp = fopen("ecg_data.txt", "r");
    if (fp == NULL) {
        printf("Error opening ecg_data.txt\n");
        return 1;
    }

    for (int i = 0; i < BUFFER_SIZE; i++) {
        if (fscanf(fp, "%f", &ecg_signal[i]) != 1) {
            printf("Error reading sample %d\n", i);
            break;
        }
    }
    fclose(fp);
    int signal_length = BUFFER_SIZE;

    // Run Pan-Tompkins pipeline
    low_pass_filter(ecg_signal, y_lp, BUFFER_SIZE);
    high_pass_filter(y_lp, y_hp, BUFFER_SIZE);
    derivative(y_hp, y_der, BUFFER_SIZE);
    square(y_der, y_sq, BUFFER_SIZE);
    integration(y_sq, y_int, BUFFER_SIZE);
    peak_detection(y_int, BUFFER_SIZE, FS, r_peaks, r_peak_binary, &num_peaks, &heart_rate);

    // Debug values (optional)
    //volatile float hr_debug = heart_rate;
    //volatile int detected_peaks = num_peaks;
    //volatile int first_peak_index = r_peaks[0];
    printf("Processing complete.\n");
    printf("Heart Rate: %.2f BPM\n", heart_rate);

    printf("R peaks detected at following indices:\n");
    for(int i = 0; i < num_peaks; i++) {
        printf("%d ", r_peaks[i]);
    }

    while (1); // Halt here for debugger and graph inspection
    return 0;
}
