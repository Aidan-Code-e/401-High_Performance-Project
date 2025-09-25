#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uthash.h"
#include <chrono>

#include <algorithm> //std::random_shuffle

#include "tracy/Tracy.hpp"

#ifndef ZoneScoped 
    #define ZoneScoped 0
#endif

typedef unsigned char byte;

#define SIGNATURE_LEN 64

int DENSITY  = 21/100;
int SIGNATURE_DENSITY = SIGNATURE_LEN * 21/100;
int PARTITION_SIZE;

int inverse[256];
const char* alphabet = "CSTPAGNDEQHRKMILVFYW";

void seed_random(char* term, int length);
short random_num(short max);
void Init();

int doc_sig[SIGNATURE_LEN];

int WORDLEN;
FILE *sig_file;

typedef struct
{
    char term[100];
    short sig[SIGNATURE_LEN];
    UT_hash_handle hh;
} hash_term;

hash_term *vocab = NULL;

short* compute_new_term_sig(char* term, short *term_sig){
    ZoneScoped;

    memset(term_sig, -1, SIGNATURE_DENSITY * sizeof(term_sig));
    memset(term_sig,  1, (SIGNATURE_DENSITY * sizeof(term_sig))/2);

    seed_random(term, SIGNATURE_LEN);
    std::random_shuffle(&term_sig[0], &term_sig[SIGNATURE_LEN], random_num);
    // printf("%d", term_sig[0]);
    return term_sig;
}

short *find_sig(char* term){
    ZoneScoped;
    hash_term *entry;
    HASH_FIND(hh, vocab, term, WORDLEN, entry);
    if (entry == NULL)    {
        entry = (hash_term*)malloc(sizeof(hash_term));
        strncpy_s(entry->term, sizeof(entry->term), term, WORDLEN);
        memset(entry->sig, 0, sizeof(entry->sig));
        compute_new_term_sig(term, entry->sig);
        HASH_ADD(hh, vocab, term, WORDLEN, entry);
    }

    return entry->sig;
}

void signature_add(char* term){
    ZoneScoped;
	short* term_sig = find_sig(term);
	for (int i=0; i<SIGNATURE_LEN; i++){
		doc_sig[i] += term_sig[i];
    }
}

int doc = 0;

void compute_signature(char* sequence, int length){
    ZoneScoped;
    memset(doc_sig, 0, sizeof(doc_sig));

    for (int i=0; i<length-WORDLEN+1; i++){
        signature_add(sequence+i);
    }

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
    const char* filename = "../small.fasta";
    // const char* filename = "../qut2.fasta";
    // const char* filename = "../qut3.fasta";
   
    WORDLEN = 3;
    PARTITION_SIZE = 16;
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
