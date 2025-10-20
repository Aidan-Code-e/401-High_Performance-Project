#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uthash.h"
#include <chrono>
#include <omp.h>

#include <iostream>
#include <cstdint>

#include <intrin.h>
#include <immintrin.h>

#include <algorithm> //std::random_shuffle

#include "tracy/Tracy.hpp"
#ifndef ZoneScoped 
    #define ZoneScoped 0
#endif

typedef unsigned char byte;

#define SIGNATURE_LEN 64

#define DENSITY 21
const int SIGNATURE_DENSITY = (SIGNATURE_LEN * DENSITY)/100;

int8_t SIGNATURE_TEMPLATE[SIGNATURE_LEN];

int PARTITION_SIZE;

int inverse[256];
const char* alphabet = "CSTPAGNDEQHRKMILVFYW";

void seed_random(char* term, int length);
short random_num(short max);
void Init();

int WORDLEN;
FILE *sig_file;

typedef struct
{
    char term[100];
    short sig[SIGNATURE_LEN];
    UT_hash_handle hh;
} hash_term;

hash_term *vocab = NULL;

void Init(){
    // This is very ugly // It's just setting the SIGNATURE_TEMPLATE array with the density of -1, 1 and 0. (it's not random though that happends later multiple times in computer_new_term_sig)
    const uint8_t step = SIGNATURE_DENSITY * sizeof(int8_t) / 2;
    memset(SIGNATURE_TEMPLATE, -1, step);
    memset(&SIGNATURE_TEMPLATE[step],  1, step);
    memset(&SIGNATURE_TEMPLATE[step*2], 0, SIGNATURE_LEN - step*2);
}

void compute_new_term_sig(char* term, int8_t *term_sig){
    ZoneScoped;

    memcpy(term_sig, SIGNATURE_TEMPLATE, SIGNATURE_LEN);
    seed_random(term, SIGNATURE_LEN);
    std::random_shuffle(&term_sig[0], &term_sig[SIGNATURE_LEN], random_num);

    return; // term_sig returned
}

/*
void compute_new_term_sig(char* term, int8_t *term_sig){
    ZoneScoped;

    const uint8_t step = SIGNATURE_DENSITY * sizeof(int8_t) / 2;

    memset(term_sig, -1, step);
    memset(&term_sig[step],  1, step);
    memset(&term_sig[step*2], 0, SIGNATURE_LEN - step*2);

    seed_random(term, SIGNATURE_LEN);
    std::random_shuffle(&term_sig[0], &term_sig[SIGNATURE_LEN], random_num);

    return; // term_sig returned
}
*/

/*
int8_t* compute_new_term_sig(char* term, int8_t *term_sig)
{
    ZoneScoped;
    memset(term_sig, 0, SIGNATURE_LEN);
    seed_random(term, WORDLEN);

    int non_zero = SIGNATURE_DENSITY;

    int positive = 0;
    while (positive < non_zero/2)
    {
        short pos = random_num(SIGNATURE_LEN);
        if (term_sig[pos] == 0){
            term_sig[pos] = 1;
            positive++;
        }
    }

    int negative = 0;
    while (negative < non_zero/2)
    {
        short pos = random_num(SIGNATURE_LEN);
        if (term_sig[pos] == 0) 
	{
            term_sig[pos] = -1;
            negative++;
        }
    }
    // printf("term: %s\n", term);
    // printf("term sig: ");
    // for(int i = 0; i<SIGNATURE_LEN; i++){
    //     printf("%d", term_sig[i]);
    // }
    // printf("\n");
    return term_sig;
}
*/

/*
int8_t *find_sig(char* term)
{
    ZoneScoped;
    hash_term *entry;
    HASH_FIND(hh, vocab, term, WORDLEN, entry);
    if (entry == NULL)
    {
        entry = (hash_term*)malloc(sizeof(hash_term));
        strncpy_s(entry->term, sizeof(entry->term), term, WORDLEN);
        memset(entry->sig, 0, sizeof(entry->sig));
        compute_new_term_sig(term, (int8_t *)entry->sig);
        HASH_ADD(hh, vocab, term, WORDLEN, entry);
    }

    return (int8_t *)entry->sig;
}
*/

// I might not be alligning right
void add_int8_simd(int8_t* a, const int8_t* b, const int8_t* c) {
    ZoneScoped;
    int i;
    for (i = 0; i + 32 <= SIGNATURE_LEN; i += 32) { // Process 32 int8_t at a time
        __m256i vec_b = _mm256_load_si256((__m256i*)&b[i]);
        __m256i vec_c = _mm256_load_si256((__m256i*)&c[i]);
        __m256i result = _mm256_add_epi8(vec_b, vec_c); // Add packed 8-bit integers
        _mm256_storeu_si256((__m256i*)&a[i], result);
    }
    // Handle remaining elements if SIGNATURE_LEN is not a multiple of 32
    for (; i < SIGNATURE_LEN; ++i) {
        a[i] = b[i] + c[i];
    }
}

int doc = 0;
void write_sig(int8_t *doc_sig){
    ZoneScoped;
    // save document number to sig file
    fwrite(&doc, sizeof(int), 1, sig_file);

    // flatten and output to sig file
    for (int i = 0; i < SIGNATURE_LEN; i += 8){
        byte c = 0;
        for (int j = 0; j < 8; j++){
            c |= (doc_sig[i+j]>0) << (7-j);
        }
        fwrite(&c, sizeof(byte), 1, sig_file);
    }
}

void compute_signature(char* sequence, int length){
    ZoneScoped;
    
    int8_t doc_sig[SIGNATURE_LEN];
    memset(doc_sig, 0, sizeof(int8_t)*SIGNATURE_LEN);

    if (length+1<WORDLEN){
        write_sig(doc_sig);
        return;
    }

    // get all signatures
    int8_t signatures[length-WORDLEN+1][SIGNATURE_LEN];

    for (int i=0; i<length-WORDLEN+1; i++){
        compute_new_term_sig(sequence+i, signatures[i]);
    }

    // add all signatures
    for (int i = 0; i<length-WORDLEN+1; i++){
        add_int8_simd(doc_sig, doc_sig, signatures[i]);
    }
    write_sig(doc_sig);
}

#define min(a,b) ((a) < (b) ? (a) : (b))

void partition(char* sequence, int length){
    ZoneScoped;
    int i=0;
    do{
        compute_signature(sequence+i, min(PARTITION_SIZE, length-i));
        i += PARTITION_SIZE/2;
    }while (i+PARTITION_SIZE/2 < length);
    doc++;
}

int power(int n, int e){
    ZoneScoped;
    int p = 1;
    for (int j=0; j<e; j++){
        p *= n;
    }
    return p;
}

int main(int argc, char* argv[])
{
    ZoneScoped;

    Init();

    const char* filename = "../small.fasta";
    // const char* filename = "../qut2.fasta";
    // const char* filename = "../qut3.fasta";

    WORDLEN = SIGNATURE_LEN;
    PARTITION_SIZE = 1024;
    int WORDS = power(20, WORDLEN);

    for (int i=0; i<strlen(alphabet); i++){
        inverse[alphabet[i]] = i;
    }
    auto start = std::chrono::high_resolution_clock::now();

    FILE* file;
    errno_t OK = fopen_s(&file, filename, "r");

    if (OK != 0){
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return 1;
    }

    char outfile[256];
    sprintf_s(outfile, 256, "%s.part%d_sigs%02d_%d", filename, PARTITION_SIZE, WORDLEN, SIGNATURE_LEN);
    fopen_s(&sig_file, outfile, "w");

    char buffer[10000];
    while (!feof(file)){
        fgets(buffer, 10000, file); // skip meta data line
        fgets(buffer, 10000, file);
        int n = (int)strlen(buffer) - 1;
        buffer[n] = 0;
        partition(buffer, n);
    }
    fclose(file);
    fclose(sig_file);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    printf("%s %f seconds\n", filename, duration.count());

    return 0;
}
