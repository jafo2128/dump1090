/****************************************************************************
*
* Name: resamp.c
*
* Synopsis: Resamples a signal.
*
* Description: See resamp.h.
*
* by Grant R. Griffin
* Copyright 2001-2015, Iowegian International Corporation
* (http://www.iowegian.com)
*
*                          The Wide Open License (WOL)
*
* Permission to use, copy, modify, distribute and sell this software and its
* documentation for any purpose is hereby granted without fee, provided that
* the above copyright notice and this license appear in all source copies. 
* THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
* ANY KIND. See http://www.dspguru.com/wide-open-license for more information.
*
*****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdint.h>

#include "fir_filter.h"
#include "fir_filter_coeffs.h"

/****************************************************************************/
void resamp1(int interp_factor_L, int decim_factor_M, int num_taps_per_phase,
             int *p_current_phase, const int16_t * const p_H,
             int8_t * const p_Z, int num_inp, const int8_t *p_inp,
             int8_t * p_out, int *p_num_out)
{
    int tap, num_out, num_new_samples, phase_num = *p_current_phase;
    const int16_t *p_coeff;
    int32_t sum;

    num_out = 0;
    while (num_inp > 0) {

        /* figure out how many new samples to shift into Z delay line */
        num_new_samples = 0;
        while (phase_num >= interp_factor_L) {
            /* decrease phase number by interpolation factor L */
            phase_num -= interp_factor_L;
            num_new_samples++;
            if (--num_inp == 0) {
                break;
            }
        }

        if (num_new_samples >= num_taps_per_phase) {
            /* the new samples are bigger than the size of Z:
               fill the entire Z with the tail of new inputs */
            p_inp += (num_new_samples - num_taps_per_phase);
            num_new_samples = num_taps_per_phase;
        }

        /* copy new samples into Z */

        /* shift Z delay line up to make room for next samples */
        for (tap = num_taps_per_phase - 1; tap >= num_new_samples; tap--) {
            p_Z[tap] = p_Z[tap - num_new_samples];
        }

        /* copy next samples from input buffer to bottom of Z */
        for (tap = num_new_samples - 1; tap >= 0; tap--) {
            p_Z[tap] = *p_inp++;
        }

        /* calculate outputs */
        while (phase_num < interp_factor_L) {
            /* point to the current polyphase filter */
            p_coeff = p_H + phase_num;

            /* calculate FIR sum */
            sum = 0.0;
            for (tap = 0; tap < num_taps_per_phase; tap++) {
                sum += *p_coeff * p_Z[tap];
                p_coeff += interp_factor_L;   /* point to next coefficient */
            }
            *p_out++ = sum >> FILTER_PRECISION;     /* store sum and point to next output */
            num_out++;

            /* decrease phase number by decimation factor M */
            phase_num += decim_factor_M;
        }
    }

    /* pass back to caller phase number (for next call) and number of
       outputs */
    *p_current_phase = phase_num;
    *p_num_out = num_out;
}

sem_t threadSync;

static struct resampThreadContext {
  int currentPhase;
  sem_t dataAvailable;
  pthread_mutex_t dataAvailableMutex;
  int8_t inBuffer[262144];
  int inLen;
  int8_t outBuffer[262144];
  int outLen;
  int8_t zBuf[FILTER_TAP_NUM];
  char *name;
} realContext, imagContext;

void *resampThread(void *arg) {
  int i;
  struct resampThreadContext *ctx = (struct resampThreadContext *)arg;

  ctx->currentPhase = 0;

  for(i = 0; i < FILTER_TAP_NUM; ctx->zBuf[i++] = 0.);

  while(1) {
    sem_wait(&ctx->dataAvailable);
    resamp1(3, 10, FILTER_TAP_NUM/3, &ctx->currentPhase,
      filter_taps, ctx->zBuf, ctx->inLen, ctx->inBuffer,
      ctx->outBuffer, &ctx->outLen
    );
    sem_post(&threadSync);
  }
  return NULL;
}

static pthread_t realThread, imagThread;
void setup_resamp_threads() {
  pthread_mutex_init(&realContext.dataAvailableMutex, NULL);
  pthread_mutex_init(&imagContext.dataAvailableMutex, NULL);
  
  pthread_mutex_lock(&realContext.dataAvailableMutex);
  pthread_mutex_lock(&imagContext.dataAvailableMutex);

  realContext.name = "real";
  imagContext.name = "imag";

  sem_init(&threadSync, 0, 0);
  sem_init(&realContext.dataAvailable, 0, 0);
  sem_init(&imagContext.dataAvailable, 0, 0);

  pthread_create(&realThread, NULL, resampThread, &realContext);
  pthread_create(&imagThread, NULL, resampThread, &imagContext);
}

#ifdef __ARM__
#include <arm_neon.h>
#else
typedef struct div {
  int32_t quot;
  int32_t rem;
} div_t;
div_t div(int32_t a, int32_t d) {
  div_t res;
  res.quot = a/d;
  res.rem = a%d;
  return res;
}
#endif

void interleave(const uint8_t *srcA, const uint8_t *srcB, uint8_t *dstAB, size_t dstABLength)
{
#if defined __ARM_NEON__
    // attempt to use NEON intrinsics

    // iterate 32-bytes at a time
    div_t dstABLength_32 = div(dstABLength, 32);
    if (dstABLength_32.rem == 0)
    {
        while (dstABLength_32.quot --> 0)
        {
            const uint8x16_t a = vld1q_u8(srcA);
            const uint8x16_t b = vld1q_u8(srcB);
            const uint8x16x2_t ab = { a, b };
            vst2q_u8(dstAB, ab);
            srcA += 16;
            srcB += 16;
            dstAB += 32;
        }
        return;
    }

    // iterate 16-bytes at a time
    div_t dstABLength_16 = div(dstABLength, 16);
    if (dstABLength_16.rem == 0)
    {
        while (dstABLength_16.quot --> 0)
        {
            const uint8x8_t a = vld1_u8(srcA);
            const uint8x8_t b = vld1_u8(srcB);
            const uint8x8x2_t ab = { a, b };
            vst2_u8(dstAB, ab);
            srcA += 8;
            srcB += 8;
            dstAB += 16;
        }
        return;
    }
#endif

    // if the bytes were not aligned properly
    // or NEON is unavailable, fall back to
    // an optimized iteration

    // iterate 8-bytes at a time
    div_t dstABLength_8 = div(dstABLength, 8);
    if (dstABLength_8.rem == 0)
    {
        typedef union
        {
            uint64_t wide;
            struct { uint8_t a1; uint8_t b1; uint8_t a2; uint8_t b2; uint8_t a3; uint8_t b3; uint8_t a4; uint8_t b4; } narrow;
        } ab8x8_t;

        uint64_t *dstAB64 = (uint64_t *)dstAB;
        int j = 0;
        for (int i = 0; i < dstABLength_8.quot; i++)
        {
            ab8x8_t cursor;
            cursor.narrow.a1 = srcA[j  ];
            cursor.narrow.b1 = srcB[j++];
            cursor.narrow.a2 = srcA[j  ];
            cursor.narrow.b2 = srcB[j++];
            cursor.narrow.a3 = srcA[j  ];
            cursor.narrow.b3 = srcB[j++];
            cursor.narrow.a4 = srcA[j  ];
            cursor.narrow.b4 = srcB[j++];
            dstAB64[i] = cursor.wide;
        }
        return;
    }

    // iterate 4-bytes at a time
    div_t dstABLength_4 = div(dstABLength, 4);
    if (dstABLength_4.rem == 0)
    {
        typedef union
        {
            uint32_t wide;
            struct { uint8_t a1; uint8_t b1; uint8_t a2; uint8_t b2; } narrow;
        } ab8x4_t;

        uint32_t *dstAB32 = (uint32_t *)dstAB;
        int j = 0;
        for (int i = 0; i < dstABLength_4.quot; i++)
        {
            ab8x4_t cursor;
            cursor.narrow.a1 = srcA[j  ];
            cursor.narrow.b1 = srcB[j++];
            cursor.narrow.a2 = srcA[j  ];
            cursor.narrow.b2 = srcB[j++];
            dstAB32[i] = cursor.wide;
        }
        return;
    }

    // iterate 2-bytes at a time
    div_t dstABLength_2 = div(dstABLength, 2);
    typedef union
    {
        uint16_t wide;
        struct { uint8_t a; uint8_t b; } narrow;
    } ab8x2_t;

    uint16_t *dstAB16 = (uint16_t *)dstAB;
    for (int i = 0; i < dstABLength_2.quot; i++)
    {
        ab8x2_t cursor;
        cursor.narrow.a = srcA[i];
        cursor.narrow.b = srcB[i];
        dstAB16[i] = cursor.wide;
    }
}

void deinterleave(const uint8_t *srcAB, uint8_t *dstA, uint8_t *dstB, size_t srcABLength)
{
#if defined __ARM_NEON__
    // attempt to use NEON intrinsics

    // iterate 32-bytes at a time
    div_t srcABLength_32 = div(srcABLength, 32);
    if (srcABLength_32.rem == 0)
    {
        while (srcABLength_32.quot --> 0)
        {
            const uint8x16x2_t ab = vld2q_u8(srcAB);
            vst1q_u8(dstA, ab.val[0]);
            vst1q_u8(dstB, ab.val[1]);
            srcAB += 32;
            dstA += 16;
            dstB += 16;
        }
        return;
    }

    // iterate 16-bytes at a time
    div_t srcABLength_16 = div(srcABLength, 16);
    if (srcABLength_16.rem == 0)
    {
        while (srcABLength_16.quot --> 0)
        {
            const uint8x8x2_t ab = vld2_u8(srcAB);
            vst1_u8(dstA, ab.val[0]);
            vst1_u8(dstB, ab.val[1]);
            srcAB += 16;
            dstA += 8;
            dstB += 8;
        }
        return;
    }
#endif

    // if the bytes were not aligned properly
    // or NEON is unavailable, fall back to
    // an optimized iteration

    // iterate 8-bytes at a time
    div_t srcABLength_8 = div(srcABLength, 8);
    if (srcABLength_8.rem == 0)
    {
        typedef union
        {
            uint64_t wide;
            struct { uint8_t a1; uint8_t b1; uint8_t a2; uint8_t b2; uint8_t a3; uint8_t b3; uint8_t a4; uint8_t b4; } narrow;
        } ab8x8_t;

        uint64_t *srcAB64 = (uint64_t *)srcAB;
        int j = 0;
        for (int i = 0; i < srcABLength_8.quot; i++)
        {
            ab8x8_t cursor;
            cursor.wide = srcAB64[i];
            dstA[j  ] = cursor.narrow.a1;
            dstB[j++] = cursor.narrow.b1;
            dstA[j  ] = cursor.narrow.a2;
            dstB[j++] = cursor.narrow.b2;
            dstA[j  ] = cursor.narrow.a3;
            dstB[j++] = cursor.narrow.b3;
            dstA[j  ] = cursor.narrow.a4;
            dstB[j++] = cursor.narrow.b4;
        }
        return;
    }

    // iterate 4-bytes at a time
    div_t srcABLength_4 = div(srcABLength, 4);
    if (srcABLength_4.rem == 0)
    {
        typedef union
        {
            uint32_t wide;
            struct { uint8_t a1; uint8_t b1; uint8_t a2; uint8_t b2; } narrow;
        } ab8x4_t;

        uint32_t *srcAB32 = (uint32_t *)srcAB;
        int j = 0;
        for (int i = 0; i < srcABLength_4.quot; i++)
        {
            ab8x4_t cursor;
            cursor.wide = srcAB32[i];
            dstA[j  ] = cursor.narrow.a1;
            dstB[j++] = cursor.narrow.b1;
            dstA[j  ] = cursor.narrow.a2;
            dstB[j++] = cursor.narrow.b2;
        }
        return;
    }

    // iterate 2-bytes at a time
    div_t srcABLength_2 = div(srcABLength, 2);
    typedef union
    {
        uint16_t wide;
        struct { uint8_t a; uint8_t b; } narrow;
    } ab8x2_t;

    uint16_t *srcAB16 = (uint16_t *)srcAB;
    for (int i = 0; i < srcABLength_2.quot; i++)
    {
        ab8x2_t cursor;
        cursor.wide = srcAB16[i];
        dstA[i] = cursor.narrow.a;
        dstB[i] = cursor.narrow.b;
    }
}

void resamp_complex(const int8_t *in, int num,
                    int8_t *out, int *num_out) {
  deinterleave((uint8_t *)in, (uint8_t *)realContext.inBuffer, (uint8_t *)imagContext.inBuffer, num);

  realContext.inLen = num/2;
  sem_post(&realContext.dataAvailable);

  imagContext.inLen = num/2;
  sem_post(&imagContext.dataAvailable);

  sem_wait(&threadSync);
  sem_wait(&threadSync);

  *num_out = realContext.outLen;
  
  interleave((uint8_t *)realContext.outBuffer, (uint8_t *)imagContext.outBuffer, (uint8_t *)out, num);

  if(realContext.outLen != imagContext.outLen) {
    fprintf(stderr, "Resampling synchronizer out of sync: %u, %u\n", realContext.outLen, imagContext.outLen);
  }
}
