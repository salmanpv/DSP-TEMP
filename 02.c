// loop variable declarartion outside the loop, not assigning and declaring variable at same time

#include <math.h>
#include <string.h>
#include <stdio.h>

#define FS 360.0f          // Sampling frequency in Hz
#define WINDOW_SIZE 30     // Integration window size
#define BUFFER_SIZE 900    // Number of ECG samples
#define MAX_PEAKS 100      // Maximum number of peaks

// Buffers
float ecg_signal[BUFFER_SIZE];
float y_lp[BUFFER_SIZE];
float y_hp[BUFFER_SIZE];
float y_der[BUFFER_SIZE];
float y_sq[BUFFER_SIZE];
float y_int[BUFFER_SIZE];
int r_peak_binary[BUFFER_SIZE];
int r_peaks[MAX_PEAKS];
float heart_rate = 0.0f;
int num_peaks = 0;

// Filters and stages
void low_pass_filter(float *x, float *y, int N) {
    int n;
    memset(y, 0, N * sizeof(float));
    for (n = 12; n < N; n++) {
        y[n] = 2 * y[n - 1] - y[n - 2] + x[n] - 2 * x[n - 6] + x[n - 12];
    }
}

void high_pass_filter(float *x, float *y, int N) {
    int n;
    memset(y, 0, N * sizeof(float));
    for (n = 32; n < N; n++) {
        y[n] = y[n - 1] - (1.0f / 32.0f) * x[n] + x[n - 16] - x[n - 17] + (1.0f / 32.0f) * x[n - 32];
    }
}

void derivative(float *x, float *y, int N) {
    int n;
    memset(y, 0, N * sizeof(float));
    for (n = 4; n < N; n++) {
        y[n] = (1.0f / 8.0f) * (2 * x[n] + x[n - 1] - x[n - 3] - 2 * x[n - 4]);
    }
}

void square(float *x, float *y, int N) {
    int n;
    for (n = 0; n < N; n++) {
        y[n] = x[n] * x[n];
    }
}

void integration(float *x, float *y, int N) {
    int n, i;
    memset(y, 0, N * sizeof(float));
    for (n = WINDOW_SIZE - 1; n < N; n++) {
        float sum = 0.0f;
        for (i = 0; i < WINDOW_SIZE; i++) {
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
    int n, i;

    memset(r_peak_binary, 0, N * sizeof(int));

    float max_val = x[0];
    for (i = 1; i < N; i++) {
        if (x[i] > max_val) max_val = x[i];
    }

    SPKI = 0.25f * max_val;
    NPKI = 0.125f * max_val;
    THRESHOLD1 = NPKI + 0.25f * (SPKI - NPKI);
    THRESHOLD2 = 0.5f * THRESHOLD1;

    for (n = 1; n < N - 1; n++) {
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
        for (i = 1; i < *num_peaks; i++) {
            float rr = (float)(r_peaks[i] - r_peaks[i - 1]) / fs;
            rr_sum += 60.0f / rr;
        }
        *heart_rate = rr_sum / (*num_peaks - 1);
    } else {
        *heart_rate = 0.0f;
    }
}

int main() {
    FILE *fp = fopen("ecg_signal.txt", "r");
    int i;
    if (fp == NULL) {
        printf("Error opening ecg_data.txt\n");
        return 1;
    }

    for (i = 0; i < BUFFER_SIZE; i++) {
        if (fscanf(fp, "%f", &ecg_signal[i]) != 1) {
            printf("Error reading sample %d\n", i);
            break;
        }
    }
    fclose(fp);

    low_pass_filter(ecg_signal, y_lp, BUFFER_SIZE);
    high_pass_filter(y_lp, y_hp, BUFFER_SIZE);
    derivative(y_hp, y_der, BUFFER_SIZE);
    square(y_der, y_sq, BUFFER_SIZE);
    integration(y_sq, y_int, BUFFER_SIZE);
    peak_detection(y_int, BUFFER_SIZE, FS, r_peaks, r_peak_binary, &num_peaks, &heart_rate);

    printf("Processing complete.\n");
    printf("Heart Rate: %.2f BPM\n", heart_rate);

    printf("R peaks detected at following indices:\n");
    for (i = 0; i < num_peaks; i++) {
        printf("%d ", r_peaks[i]);
    }

    while (1); // Stay here for debugging
    return 0;
}
