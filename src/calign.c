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

void destroy_seq_pair(seq_pair_t pair)
{
    free(pair->a);
    free(pair->b);

    free(pair);

    return;
}

int local_alignment_score(seq_pair_t problem)
{
    unsigned int n = problem->alen+1;
    unsigned int m = problem->blen+1;
    int score_arr[m];
    int curr_val;
    unsigned int i, j;
    int last_sub_val;
    int curr_sub_val;
    int curr_ins_val;
    int curr_del_val;

    int max_score = 0;
    int result_row_index = 0;
    int result_col_index = 0;

    for (i = 0; i<m; i++) {
        score_arr[i] = 0;
    }

    for (i = 1; i<n; i++) {
        last_sub_val = 0;
        for (j = 1; j<m; j++) {
            int nw_score = (strncmp(problem->a+i-1, problem->b+j-1, 1)==0) ? MATCH : MISMATCH;

            curr_val = INT_MIN;

            curr_del_val = score_arr[j]+ GAP;
            curr_ins_val = score_arr[j-1]+ GAP;
            curr_sub_val = last_sub_val+nw_score;
            curr_val = (curr_ins_val<curr_del_val) ? curr_del_val :
                       curr_ins_val;
            curr_val = (curr_val<curr_sub_val) ? curr_sub_val : curr_val;
            curr_val = (curr_val<0) ? 0 : curr_val;

            if (max_score<curr_val) {
                max_score = curr_val;
                result_row_index = i;
                result_col_index = j;
            }
            last_sub_val = score_arr[j];
            score_arr[j] = curr_val;
        }
    }
    return max_score;
}

int semiglobal_alignment_score(seq_pair_t problem, int top, int left, int right, int bottom)
{
    unsigned int n = problem->alen+1;
    unsigned int m = problem->blen+1;
    int score_arr[m];
    int curr_val;
    unsigned int i, j;
    int last_sub_val;
    int curr_sub_val;
    int curr_ins_val;
    int curr_del_val;

    int max_score = 0;
    int result_row_index = 0;
    int result_col_index = 0;

    for (i = 0; i<m; ++i) {
        score_arr[i] = left ? 0 : i*GAP; //non free begin gaps of vertical sequence
    }

    for (i = 1; i<n; ++i) {
        last_sub_val = top ? 0 : (i-1)*GAP; // (non) free begin gaps of horizontal sequence
        for (j = 1; j<m; ++j) {
            int nw_score = (strncmp(problem->a+i-1, problem->b+j-1, 1)==0) ? MATCH : MISMATCH;

            curr_val = INT_MIN;

            curr_del_val = score_arr[j]+ GAP;
            curr_ins_val = score_arr[j-1]+ GAP;
            curr_sub_val = last_sub_val+nw_score;
            curr_val = (curr_ins_val<curr_del_val) ? curr_del_val :
                       curr_ins_val;
            curr_val = (curr_val<curr_sub_val) ? curr_sub_val : curr_val;

            last_sub_val = score_arr[j];
            score_arr[j] = curr_val;
        }
        if (bottom==1) {
            if (max_score<score_arr[m-1]) {
                max_score = score_arr[m-1];
                result_row_index = i;
                result_col_index = m-1;
            }
        }
    }
    if (right==1) {
        for (j = 0; j<m; ++j) {
            if (max_score<score_arr[j]) {
                max_score = score_arr[j];
                result_row_index = n-1;
                result_col_index = j;
            }
        }
    }

    return max_score;
}

int banded_semiglobal_alignment_score(seq_pair_t problem, int top, int left, int right, int bottom, int lower_diag,
        int upper_diag)
{
    unsigned int m = problem->alen+1;
    unsigned int n = problem->blen+1;
    int array_size = MIN(n, upper_diag-lower_diag+1);

    if (n==array_size) return semiglobal_alignment_score(problem, top, left, right, bottom);

    int score_arr[array_size];
    int curr_val;
    unsigned int i, j;

    int curr_sub_val;
    int curr_ins_val;
    int curr_del_val;
    int sub_score;
    int k;
    int start_index;
    int a_str_index;
    int max_score = 0;
    int result_row_index = 0;
    int result_col_index = 0;

    for (i = 0; i<array_size; ++i) {
        score_arr[i] = top ? 0 : i*GAP; //non free begin gaps of horizontal sequence
    }

    for (i = 1; i<m; ++i) {
        k = 0;
        start_index = (1-lower_diag-i);
        a_str_index = -start_index;
        if (0>=start_index) start_index = 0;
        if (a_str_index<0) a_str_index = 0;

        for (j = start_index; j<MIN(array_size, n-i-lower_diag); ++j) {

            //set substitution value
//            sub_score = (strncmp(problem->b+i-1, problem->a+a_str_index, 1)==0) ? MATCH : MISMATCH;
            sub_score = (*(problem->b+i-1)==*(problem->a+a_str_index)) ? MATCH : MISMATCH;
            curr_sub_val = score_arr[j]+sub_score;

            //set insertion value
            if (k>0) {
                curr_ins_val = score_arr[j-1]+ GAP;
            }
            else if (k==0) {
                if (i<=-lower_diag) curr_ins_val = left ? 0 : i*GAP;
                else curr_ins_val = INT_MIN;
            }

            //set deletion value
            if (j==array_size-1) curr_del_val = INT_MIN;
            else curr_del_val = score_arr[j+1]+ GAP;

            curr_val = MAX(curr_sub_val, MAX(curr_del_val, curr_ins_val));

            score_arr[j] = curr_val;
            k++;
            a_str_index++;
        }

        if (right==1 && i+upper_diag+1>=n) {
            if (max_score<score_arr[array_size-1]) {
                max_score = score_arr[array_size-1];
                result_row_index = i;
                result_col_index = array_size-1;
            }
        }
    }
    if (bottom==1) {
        for (j = 0; j<array_size; ++j) {
            if (max_score<score_arr[j]) {
                max_score = score_arr[j];
                result_row_index = m-1;
                result_col_index = j;
            }
        }
    }

    return max_score;
}

int alignment_score(seq_pair_t problem, char* alignment_type)
{
    if (strcmp(alignment_type, "local")==0) return local_alignment_score(problem);

    if (strcmp(alignment_type, "semiglobal")==0) return semiglobal_alignment_score(problem, TOP, LEFT, RIGHT, BOTTOM);
}

int alignment_score2(const char* a, int alen, const char* b, int blen)
{
    seq_pair problem;
    problem.a = a;
    problem.alen = alen;
    problem.b = b;
    problem.blen = blen;
    return semiglobal_alignment_score(&problem, 1, 0, 1, 1);
}

void alignment_score_all(const char* filename, int len, const char** strings, char* alignment_type)
{
    FILE* fp;
    fp = fopen(filename, "w");
    fprintf(fp, "%i\n", len);
    int max_score;
    seq_pair problem;
    for (int i = 0; i<len; ++i) {
        for (int j = 0; j<len; ++j) {
            if (i!=j) {
                problem.a = strings[i];
                problem.alen = strlen(problem.a);
                problem.b = strings[j];
                problem.blen = strlen(problem.b);
                max_score = alignment_score(&problem, alignment_type);
                if (max_score>0) fprintf(fp, "%i\t%i\t%i\n", i, j, max_score);
            }
        }
    }
    fclose(fp);
}

void banded_alignment_score_all(const char* filename, int len, const char** strings, int lower_diag, int upper_diag)
{
    FILE* fp;
    fp = fopen(filename, "w");
    fprintf(fp, "%i\n", len);
    int max_score;
    seq_pair problem;
    for (int i = 0; i<len; ++i) {
        for (int j = 0; j<len; ++j) {
            if (i!=j) {
                problem.a = strings[i];
                problem.alen = strlen(problem.a);
                problem.b = strings[j];
                problem.blen = strlen(problem.b);
                max_score = banded_semiglobal_alignment_score(&problem, TOP, LEFT, RIGHT, BOTTOM, lower_diag,
                        upper_diag);
                if (max_score>0) fprintf(fp, "%i\t%i\t%i\n", i, j, max_score);
            }
        }
    }
    fclose(fp);
}

int main(int argc, const char** argv)
{

//    if (argc!=3) {
//        printf("usage: calign SEQ1 SEQ2\n");
//        exit(1);
//    }

    {
        char b[] = "TCTTGCAACATGCTTATGTAACGATGAGTTAGCAATATGCCTTACAAGGAAAGAAAAGGCACCGTGCATGCCGATTGGTGGTAGTAAGGTGGTACGATCG";
        char c[] = "TTATGTAACGATGAGTTAGCAATATGCCTTACAAGGAAAGAAAAGGCACCGTGCATGCCGATTGGTGGTAGTAAGGTGGTACGATCGTGCCTTATTAGGA";
//        char c[] = "GTGGGTATCAGATATCAGA";
//        char b[] =       "ATCAGATATCAGAAAAAAA";
        seq_pair foo;
        foo.a = b;
        foo.b = c;
        foo.alen = strlen(b);
        foo.blen = strlen(c);
        printf("%d\n", alignment_score(&foo, "semiglobal"));
        printf("%d", banded_semiglobal_alignment_score(&foo, TOP, LEFT, RIGHT, BOTTOM, 0, 25));


    }

    exit(0);
}
