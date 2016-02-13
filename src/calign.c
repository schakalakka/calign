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

#include "calign.h"


/* reverse a string in place, return str */
static char *reverse(char *str) {
    char *left = str;
    char *right = left + strlen(str) - 1;
    char tmp;

    while (left < right) {
        tmp = *left;
        *(left++) = *right;
        *(right--) = tmp;
    }

    return str;
}

// works globally
static seq_pair_t traceback(seq_pair_t problem, matrix_t S, bool local) {
    seq_pair_t result = malloc(sizeof(seq_pair));
    unsigned int i = S->m - 1;
    unsigned int j = S->n - 1;
    unsigned int k = 0;
    char c[S->m + S->n + 1];
    char d[S->m + S->n + 1];

    memset(c, '\0', sizeof(c));
    memset(d, '\0', sizeof(d));

    if (local == true) {
        unsigned int l, m;
        double max = FLT_MIN;

        for (l = 0; l < S->m; l++) { //TODO S->m? shouldn't it be n?
            for (m = 0; m < S->n; m++) {
                if (S->mat[l][m]->score > max) {
                    i = l;
                    j = m;
                }
            }
        }
    }

    if (S->mat[i][j]->prev[0] != 0 && S->mat[i][j]->prev[1] != 0) {
        while (i > 0 || j > 0) {
            unsigned int new_i = S->mat[i][j]->prev[0];
            unsigned int new_j = S->mat[i][j]->prev[1];

            if (new_i < i)
                *(c + k) = *(problem->a + i - 1);
            else
                *(c + k) = '-';

            if (new_j < j)
                *(d + k) = *(problem->b + j - 1);
            else
                *(d + k) = '-';

            k++;

            i = new_i;
            j = new_j;
        }
    }

    result->a = malloc(sizeof(char) * k + 1);
    result->b = malloc(sizeof(char) * k + 1);

    memset(result->a, '\0', sizeof(result->a));
    memset(result->b, '\0', sizeof(result->b));

    reverse(c);
    reverse(d);

    strcpy(result->a, c);
    strcpy(result->b, d);

    result->alen = k;
    result->blen = k;

    return result;
}

static matrix_t create_matrix(unsigned int m, unsigned int n) {
    matrix_t S = malloc(sizeof(matrix));
    unsigned int i, j;

    S->m = m;
    S->n = n;

    S->mat = malloc(sizeof(entry_t) * m * n);

    for (i = 0; i < m; i++) {
        S->mat[i] = malloc(sizeof(entry_t) * n);
    }

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            S->mat[i][j] = malloc(sizeof(entry));
        }
    }

    return S;
}

void destroy_matrix(matrix_t S) {
    unsigned int i, j;

    for (i = 0; i < S->m; i++) {
        for (j = 0; j < S->n; j++) {
            free(S->mat[i][j]);
        }
    }

    free(S);

    return;
}

void destroy_seq_pair(seq_pair_t pair) {
    free(pair->a);
    free(pair->b);

    free(pair);

    return;
}

static seq_pair_t smith_waterman(seq_pair_t problem, bool local) {
    unsigned int m = problem->alen + 1;
    unsigned int n = problem->blen + 1;
    matrix_t S = create_matrix(m, n);
    seq_pair_t result;
    unsigned int i, j, k, l;

    S->mat[0][0]->score = 0;
    S->mat[0][0]->prev[0] = 0;
    S->mat[0][0]->prev[1] = 0;

    for (i = 1; i <= problem->alen; i++) {
        S->mat[i][0]->score = 0;
        S->mat[i][0]->prev[0] = 0;//i - 1;
        S->mat[i][0]->prev[1] = 0;
    }

    for (j = 1; j <= problem->blen; j++) {
        S->mat[0][j]->score = 0;
        S->mat[0][j]->prev[0] = 0;
        S->mat[0][j]->prev[1] = 0;//j - 1;
    }

    for (i = 1; i <= problem->alen; i++) {
        for (j = 1; j <= problem->blen; j++) {
            int nw_score = (strncmp(problem->a + (i - 1), problem->b + (j - 1), 1) == 0) ? MATCH : MISMATCH;

            S->mat[i][j]->score = INT_MIN;
            S->mat[i][j]->prev[0] = 0;
            S->mat[i][j]->prev[1] = 0;

            for (k = 0; k <= 1; k++) {
                for (l = 0; l <= 1; l++) {
                    int val;

                    if (k == 0 && l == 0) {
                        continue;
                    } else if (k > 0 && l > 0) {
                        val = nw_score;
                    } else if (k > 0 || l > 0) {
                        val = GAP;
                    } else {
                        // do nothing..
                    }

                    val += S->mat[i - k][j - l]->score;

                    if (val > S->mat[i][j]->score) {
                        S->mat[i][j]->score = val;
                        S->mat[i][j]->prev[0] = i - k;
                        S->mat[i][j]->prev[1] = j - l;
                    }
                }
            }
        }
    }

    printf("%i\t%i\n", m, n);
    int x, y;
    for (x = 0; x < m; ++x) {
        for (y = 0; y < n; y++) {
            printf("%d\t", S->mat[x][y]->score);
        }
        printf("\n");
    }

    for (x = 0; x < m; ++x) {
        for (y = 0; y < n; y++) {
            printf("(%d,%d)\t", S->mat[x][y]->prev[0], S->mat[x][y]->prev[1]);
        }
        printf("\n");
    }

    result = traceback(problem, S, local);

    destroy_matrix(S);

    return result;
}

int *local_alignment(seq_pair_t problem) {
    seq_pair_t result;
    unsigned int n = problem->alen + 1;
    unsigned int m = problem->blen + 1;
    int *max_score = malloc(sizeof(int) * 3);
    int score_arr[m];
    int curr_val;
    unsigned int i, j;
    int last_sub_val;
    int curr_sub_val;
    int curr_ins_val;
    int curr_del_val;

    max_score[0] = 0;
    max_score[1] = 0;
    max_score[2] = 0;

    for (i = 0; i < m; i++) {
        score_arr[i] = 0;
    }


    for (i = 1; i < n; i++) {
        last_sub_val = 0;
        for (j = 1; j < m; j++) {
            int nw_score = (strncmp(problem->a + i - 1, problem->b + j - 1, 1) == 0) ? MATCH : MISMATCH;

            curr_val = INT_MIN;

            curr_del_val = score_arr[j] + GAP;
            curr_ins_val = score_arr[j - 1] + GAP;
            curr_sub_val = last_sub_val + nw_score;
            curr_val = (curr_ins_val < curr_del_val) ? curr_del_val :
                       curr_ins_val;
            curr_val = (curr_val < curr_sub_val) ? curr_sub_val : curr_val;
            curr_val = (curr_val < 0) ? 0 : curr_val;

            if (max_score[0] < curr_val) {
                max_score[0] = curr_val;
                max_score[1] = i;
                max_score[2] = j;
            }
            last_sub_val = score_arr[j];
            score_arr[j] = curr_val;
        }
    }
    return max_score;
}

int *semiglobal_alignment(seq_pair_t problem) {
    seq_pair_t result;
    unsigned int n = problem->alen + 1;
    unsigned int m = problem->blen + 1;
    int *max_score = malloc(sizeof(int) * 3);
    int score_arr[m];
    int curr_val;
    unsigned int i, j;
    int last_sub_val;
    int curr_sub_val;
    int curr_ins_val;
    int curr_del_val;

    max_score[0] = 0;
    max_score[1] = 0;
    max_score[2] = 0;

    for (i = 0; i < m; i++) {
        score_arr[i] = 0;
    }

    for (i = 1; i < n; i++) {
        last_sub_val = 0;
        for (j = 1; j < m; j++) {
            int nw_score = (strncmp(problem->a + i - 1, problem->b + j - 1, 1) == 0) ? MATCH : MISMATCH;

            curr_val = INT_MIN;

            curr_del_val = score_arr[j] + GAP;
            curr_ins_val = score_arr[j - 1] + GAP;
            curr_sub_val = last_sub_val + nw_score;
            curr_val = (curr_ins_val < curr_del_val) ? curr_del_val :
                       curr_ins_val;
            curr_val = (curr_val < curr_sub_val) ? curr_sub_val : curr_val;

            last_sub_val = score_arr[j];
            score_arr[j] = curr_val;
        }
        if (max_score[0] < score_arr[m - 1]) {
            max_score[0] = score_arr[m - 1];
            max_score[1] = i;
            max_score[2] = m - 1;
        }
    }
    for (j = 0; j < m; j++) {
        if (max_score[0] < score_arr[j]) {
            max_score[0] = score_arr[j];
            max_score[1] = n - 1;
            max_score[2] = j;
        }
    }
    return max_score;
}

int *global_alignment(seq_pair_t problem) {
    seq_pair_t result;
    unsigned int n = problem->alen + 1;
    unsigned int m = problem->blen + 1;
    int *max_score = malloc(sizeof(int) * 3);
    int score_arr[m];
    int curr_val;
    unsigned int i, j;
    int last_sub_val;
    int curr_sub_val;
    int curr_ins_val;
    int curr_del_val;

    max_score[0] = 0;
    max_score[1] = 0;
    max_score[2] = 0;

    for (i = 0; i < m; i++) {
        score_arr[i] = i * GAP;
    }

    for (i = 1; i < n; i++) {
        last_sub_val = (i - 1) * GAP;
        score_arr[0] = i * GAP;
        for (j = 1; j < m; j++) {
            int nw_score = (strncmp(problem->a + i - 1, problem->b + j - 1, 1) == 0) ? MATCH : MISMATCH;

            curr_val = INT_MIN;

            curr_del_val = score_arr[j] + GAP;
            curr_ins_val = score_arr[j - 1] + GAP;
            curr_sub_val = last_sub_val + nw_score;
            curr_val = (curr_ins_val < curr_del_val) ? curr_del_val :
                       curr_ins_val;
            curr_val = (curr_val < curr_sub_val) ? curr_sub_val : curr_val;

            last_sub_val = score_arr[j];
            score_arr[j] = curr_val;
        }
    }
    max_score[0] = score_arr[m - 1];
    max_score[1] = n - 1;
    max_score[2] = m - 1;

    return max_score;
}


int *alignment(seq_pair_t problem, char *foo) {
    if (strcmp(foo, "local") == 0) return local_alignment(problem);

    if (strcmp(foo, "semiglobal") == 0) return semiglobal_alignment(problem);

    if (strcmp(foo, "global") == 0) return global_alignment(problem);

}


int main(int argc, const char **argv) {

    if (argc != 3) {
        printf("usage: calign SEQ1 SEQ2\n");
        exit(1);
    }

    {
        seq_pair problem;
        seq_pair_t result;
        int *max_score;
        char c[strlen(argv[1])], d[strlen(argv[2])];

        strcpy(c, argv[1]);
        strcpy(d, argv[2]);


        problem.a = c;
        problem.alen = strlen(problem.a);
        problem.b = d;
        problem.blen = strlen(problem.b);

        //result = smith_waterman(&problem, false);
        //printf("%i\n%i\n", problem.alen, problem.blen);

        max_score = alignment(&problem, "global");
        printf("Global Alignment:\t %i\t%i\t%i\n", max_score[0], max_score[1], max_score[2]);
        max_score = alignment(&problem, "semiglobal");
        printf("Semiglobal Alignment:\t %i\t%i\t%i\n", max_score[0], max_score[1], max_score[2]);
        max_score = alignment(&problem, "local");
        printf("Local Alignment:\t %i\t%i\t%i\n", max_score[0], max_score[1], max_score[2]);
    }

    exit(0);
}