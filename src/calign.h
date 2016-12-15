/*
 *  Copyright (c) 2010 Nicolaus Lance Hepler
 *
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use,
 *  copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following
 *  conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *  OTHER DEALINGS IN THE SOFTWARE.
 */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <sys/param.h>

#define MIN(X, Y) (((X) < (Y) ? (X) : (Y)))
#define MAX(X, Y) (((X) > (Y) ? (X) : (Y)))

#define GAP -5000
#define MATCH 1
#define MISMATCH -5000

#define BOOL int
#define FALSE 0
#define TRUE 1

#define TOP TRUE
#define LEFT FALSE
#define RIGHT TRUE
#define BOTTOM TRUE


typedef struct {
    char *a;
    unsigned int alen;
    char *b;
    unsigned int blen;
} seq_pair;
typedef seq_pair *seq_pair_t;


void destroy_seq_pair(seq_pair_t pair);

static seq_pair_t alignment(seq_pair_t problem, char *alignment_type);

int local_alignment_score(seq_pair_t problem);

int semiglobal_alignment_score(seq_pair_t problem, int, int, int, int);

int banded_semiglobal_alignment_score(seq_pair_t problem, int, int, int, int, int, int);

int alignment_score(seq_pair_t problem, char* alignment_type);

int alignment_score2(const char*, int, const char*, int);

void alignment_score_all(const char*, int, const char**, char*);

void banded_alignment_score_all(const char*, int len, const char** strings, int, int);




