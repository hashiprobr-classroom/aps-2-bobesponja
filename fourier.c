#include <math.h>

#include "fourier.h"

void normalize(complex s[], int n) {
    for (int k = 0; k < n; k++) {
        s[k].a /= n;
        s[k].b /= n;
    }
}

void nft(complex s[], complex t[], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k].a = 0;
        t[k].b = 0;

        for (int j = 0; j < n; j++) {
            double x = sign * 2 * PI * k * j / n;

            double cosx = cos(x);
            double sinx = sin(x);

            t[k].a += s[j].a * cosx - s[j].b * sinx;
            t[k].b += s[j].a * sinx + s[j].b * cosx;
        }
    }
}

void nft_forward(complex s[], complex t[], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(complex t[], complex s[], int n) {
    nft(t, s, n, 1);
    normalize(s, n);
}
void banana(complex s[], complex t[], int n, int sign){
    if(n == 0){
        return;
    }

    t[n-1].a = 0;
    t[n-1].b = 0;
    
    for (int j = 0; j < n; j++) {
        double x = sign * 2 * PI * (n-1) * j / n;

        double cosx = cos(x);
        double sinx = sin(x);

        t[n-1].a += s[j].a * cosx - s[j].b * sinx;
        t[n-1].b += s[j].a * sinx + s[j].b * cosx;
    }
    banana(s, t, n-1, sign);
}
void fft(complex s[], complex t[], int n, int sign) {
        complex lista_s_1[n/2];
        complex lista_s_2[n/2];
        complex lista_t_1[n/2];
        complex lista_t_2[n/2];
        int cont = 0;
        int cont_1 = 1;
            for(int i = 0; i<(n/2); i++){
            lista_s_1[i] = s[cont];
            cont+=2;
        }
        for(int j = 0; j<(n/2); j++){
            lista_s_2[j] = s[cont_1];
            cont_1+=2;
        }
        banana(lista_s_1, lista_t_1, (n/2), sign);
        banana(lista_s_2, lista_t_2, (n/2), sign);

        for (int k = 0; k < n / 2; k++) {
            double x = sign * 2 * PI * k / n;
            double cosx = cos(x);
            double sinx = sin(x);

            double wk_ti_a = lista_t_2[k].a * cosx - lista_t_2[k].b * sinx;
            double wk_ti_b = lista_t_2[k].a * sinx + lista_t_2[k].b * cosx;

            t[k].a = lista_t_1[k].a + wk_ti_a;
            t[k].b = lista_t_1[k].b + wk_ti_b;

            t[k + n / 2].a = lista_t_1[k].a - wk_ti_a;
            t[k + n / 2].b = lista_t_1[k].b - wk_ti_b;
        }
    // usar duas listas uma so com os pares e outra so com impares(indices) da lista original
    //  usar o negocio do cosseno e seno que nem no nft pq n podemos usar exp
}

void fft_forward(complex s[], complex t[], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(complex t[], complex s[], int n) {
    fft(t, s, n, 1);
    normalize(s, n);
}

void fft_forward_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
}

void filter(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;
            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x].a = g * input[y][x].a;
            output[y][x].b = g * input[y][x].b;
        }
    }
}

void filter_lp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
