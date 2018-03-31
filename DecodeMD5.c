#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include "md5.c"

void divide_function(int * dividend, int dividend_l, int * divisor, int divisor_l, int * result, int * result_l);
void multiply_function_val(int * factor1, int factor1_l, int val, int * result, int * result_l);
void multiply_function(int * factor1, int factor1_l, int * factor2, int factor2_l, int * result, int * result_l);
void add_function(int * addend1, int addend1_l, int * addend2, int addend2_l, int * result, int * result_l);
int b_base_to_decimal(int * arr, int arr_l);
void decimal_to_b_base(int val, int *result, int *result_l);
int less_than(int *arr1, int arr1_l, int * arr2, int arr2_l);
int  find_plain_text(int rank, int num_p, char * hashed_data, int *start, int s_l, int * end, int e_l, char * result, char * plain_text, int *from);
void compute_md5(char *str, unsigned char digest[16]);
void add_zeros_to_front(int * arr, int arr_l, int n);
int determine_mode(char * mode);
void arr_to_string(int * arr, int arr_l, char * result);
void usage(char * s);
void rank0(int rank, int num_p, int* first_bound, char* hash_string, int n);
void ranki(int rank, int num_p);
int init(char argc, char** args, int num_p, int* first_bound, char* hash_string, int *n, char * message);

int base;
char mapping[62];

int main(int argc, char **args)
{
    MPI_Init(NULL, NULL);

    int rank, num_p, total_time = 0, time, iteration = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_p);

    int first_bound[50], n, flag;
    char hash_string[50], message[500];
    for (int i = 0; i < iteration; i++){
        MPI_Barrier(MPI_COMM_WORLD);
        time = -1 * MPI_Wtime();

        if(rank == 0){
            flag = init(argc, args, num_p, first_bound, hash_string, &n, message);
            MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(flag){
                printf("%s", message);
            }
            else{
                rank0(rank, num_p, first_bound, hash_string, n);
            }
        }
        else{
            MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(flag == 0){
                ranki(rank, num_p);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        time += MPI_Wtime();
        total_time += time;
        if(rank == 0){
        printf("Time: %d seconds\n", time);
    }
    }
    if(rank == 0){
        printf("Average time: %.2f seconds\n", 1.0 * total_time/iteration);
    }

    MPI_Finalize();
}

int init(char argc, char** args, int num_p, int* first_bound, char* hash_string, int *n, char* message){
    int flag, i, len;
    char mode[10], tmp[10];


    if(argc != 4){
        char usg[500];
        usage(usg);
        sprintf(message, "%s", usg);
        return -1;
    }

    strcpy(mode, args[1]);
    strcpy(tmp, args[2]);
    flag = 0;
    len = strlen(tmp);
    for(i = 0; i < len; i++){
        if(tmp[i] < '0' || tmp[i] > '9')
            flag = 1;
    }
    if(flag){
        usage(message);
        return -1;
    }
    else
        *n = atoi(tmp);
    strcpy(hash_string, args[3]);
    if(strlen(hash_string) != 32){
        sprintf(message, "Invalid hash string\n");
        return -1; 
    }

    flag = determine_mode(mode);
    if(flag){
        sprintf(message, "Invalid mode\n");
        return -1;
    }

    int max[50];
    for(i = 0; i < *n; i++){
        max[i] = base - 1;
    }

    int fb_l, divisor[1];
    divisor[0] = num_p;
    divide_function(max, *n, divisor, 1, first_bound, &fb_l);
    add_zeros_to_front(first_bound, fb_l, *n);

    return 0;
}

void rank0(int rank, int num_p, int* first_bound, char* hash_string, int n){
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&base, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(first_bound, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(mapping, base, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(hash_string, 33, MPI_CHAR, 0, MPI_COMM_WORLD);

    int start[50];
    add_zeros_to_front(start, 0, n);
    char result[60], plain_text[50];
    int from, flag;
    flag = find_plain_text(rank, num_p, hash_string, start, n, first_bound, n, result, plain_text, &from);
    if(flag == 1)
        printf("Found: %s by rank %d\n", plain_text, from);
    else{
        printf("Not found...\n");
    }
}

void ranki(int rank, int num_p){
    int first_bound[50], n;
    char hash_string[50], result[50];

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&base, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(first_bound, n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(mapping, base, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(hash_string, 33, MPI_CHAR, 0, MPI_COMM_WORLD);

    int start[50], end[50], s_l, e_l, factor[50], ft_l;
    decimal_to_b_base(rank, factor, &ft_l);
    multiply_function(first_bound, n, factor, ft_l, start, &s_l);
    add_zeros_to_front(start, s_l, n);
    if(rank == num_p - 1){
        end[0] = 1;
        for(int i = 1; i < n +1; i++){
            end[i] = 0;
        }
        find_plain_text(rank, num_p, hash_string, start, n, end, n + 1, result, NULL, NULL);
    }
    else{
        int add_unit[1] = {1};
        add_function(factor, ft_l, add_unit, 1, factor, &ft_l);
        multiply_function(first_bound, n, factor, ft_l, end, &e_l);
        add_zeros_to_front(end, e_l, n);
        find_plain_text(rank, num_p, hash_string, start, n, end, n, result, NULL, NULL);
    }
}

void usage(char * s){
    sprintf(s, "Usage: <mode> <n> <hash string>\n"
    "   mode: a string combined of \'0\' for the use of digits, \'a\' for the use of normal,\n" 
    "and \'A\' for the use of capitals. \n"
    "   n is the number of characters in the plaintext. \n"
    "   hash string is md5 hashed data.\n");
}

void arr_to_string(int * arr, int arr_l, char * result){
    int i;
    for(i = 0; i < arr_l; i++){
        result[i] = mapping[arr[i]];
    }
    result[i] = '\0';
}


int determine_mode(char * mode){
    int len = strlen(mode), marks[3] = {0, 0, 0};
    if(len > 3){
        return 1;
    }
    for (int i = 0; i < len; i++){
        if(mode[i] == '0'){
            marks[0] = 1;
        }
        else if(mode[i] == 'a'){
            marks[1] = 1;
        }
        else if(mode[i] == 'A'){
            marks[2] = 1;
        }
        else{
            return 1;
        }
    }
    int i = 0, n = 0;
    if(marks[0] == 1){
        n += 10;
        for(char c = '0'; c <= '9'; c++){
            mapping[i] = c;
            i++;
        }
    }
    if(marks[1] == 1){
        n += 26;
        for(char c = 'a'; c <= 'z'; c++){
            mapping[i] = c;
            i++;
        }
    }
    if(marks[2] == 1){
        n += 26;
        for(char c = 'A'; c <= 'Z'; c++){
            mapping[i] = c;
            i++;
        }
    }
    mapping[i] = '\0';
    base = n;
    return 0;
}

void add_zeros_to_front(int * arr, int arr_l, int n){
    for(int i = 0; i < arr_l; i++){
        arr[n - 1 - i] = arr[arr_l - 1 - i];
    }
    for(int i = 0; i < n - arr_l; i++){
        arr[i] = 0;
    }
}

int find_plain_text(int rank, int num_p, char * hashed_data, int *start, int s_l, int * end, int e_l, char * result, char* found_text, int* from){//range from start to end -1
    char host[50], str_start[50], str_end[50];
    int tmp_len;
    arr_to_string(start, s_l, str_start);
    arr_to_string(end, e_l, str_end);
    MPI_Get_processor_name(host, &tmp_len);
    printf("Rank %d in %s: finding from %s to %s\n", rank, host,str_start, str_end);
    
    int prediction[50], n = s_l; 
    for(int i = 0; i < n; i++){
        prediction[i] = start[i];
    }
    char hash_output[33], text[50], plain_text[50];
    int flag;
    MPI_Request mpi_request;
    MPI_Status mpi_status;
    MPI_Irecv(plain_text, 50, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpi_request);
    while(less_than(prediction, n, end, e_l)){
        MPI_Test(&mpi_request, &flag, &mpi_status);
        if(flag){
            if(rank == 0){
                strcpy(found_text, plain_text);
                *from = mpi_status.MPI_SOURCE;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }
        arr_to_string(prediction, n, text);
        //printf("rank %d: %s\n", rank, text);
        md5_string(text, hash_output);
        if (strcmp(hashed_data, hash_output) == 0){
            if(rank == 0){
                strcpy(found_text, text);
                *from = 0;
            }
            strcpy(result, text);
            for(int i = 0; i < num_p; i++)
                if(i != rank){
                    MPI_Isend(result, n, MPI_CHAR, i, 0, MPI_COMM_WORLD, &mpi_request);
                    MPI_Wait(&mpi_request, &mpi_status);
                }
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }
        int add_unit[1] = {1}, _;
        add_function(prediction, n, add_unit, 1, prediction, &n);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Test(&mpi_request, &flag, &mpi_status);
    char x[1] = "";
    if(rank == 0){
        if(flag){
            strcpy(found_text, plain_text);
            *from = mpi_status.MPI_SOURCE;
            return 1;
        }
    }
    return 0;
}

void decimal_to_b_base(int val, int *result, int *result_l){
    int reverse_result[50], i = 0;
    while(val > 0){
        reverse_result[i] = val % base;
        val /= base;
        i++;
    }
    for(int j = 0; j < i; j++){
        result[j] = reverse_result[i - 1 - j];
    }
    *result_l = i;
}

void divide_function(int * dividend, int dividend_l, int * divisor, int divisor_l, int * result, int * result_l){
    int tmp, i_dividend, i_result = 0, divisor_val = b_base_to_decimal(divisor, divisor_l);
    if (less_than(dividend, dividend_l, divisor, divisor_l)){
        result[0] = 0;
        *result_l = 1;
        return;
    }
    i_dividend = divisor_l - 1;
    *result_l = 0;
    tmp = b_base_to_decimal(dividend, divisor_l);
    while(i_dividend < dividend_l){
        result[i_result] = tmp / divisor_val;
        i_result++;
        tmp = tmp % divisor_val;
        i_dividend++;
        tmp = tmp * base + dividend[i_dividend];
    }
    *result_l = i_result;
}

int less_than(int *arr1, int arr1_l, int * arr2, int arr2_l){
    if(arr1_l < arr2_l)
        return 1;
    if(arr1_l > arr2_l)
        return 0;
    for(int i = 0; i < arr1_l; i++){
        if(arr1[i] < arr2[i]){
            return 1;
        }
        else if(arr1[i] > arr2[i]){
            return 0;
        }
    }
    return 0;
}

int b_base_to_decimal(int *arr, int arr_l){
    int result = 0, e = arr_l - 1; 
    for(int i = 0; i < arr_l; i++){
        result += arr[i] * (int)pow(base, e);
        e--;
    }
    return result;
}

void multiply_function_val(int * factor, int factor_l, int val, int * result, int * result_l){
    int r = 0, tmp_result[factor_l + 1];
    int i_tmp_result = 0, t;
    for(int i = factor_l - 1; i >= 0; i--){
        t = val * factor[i] + r;
        tmp_result[i_tmp_result] = t % base;
        r = t / base;
        i_tmp_result++;
    }
    for(int i = 0; i < i_tmp_result; i++){
        result[i] = tmp_result[i_tmp_result - i - 1];
    }
    *result_l = i_tmp_result;
}

void multiply_function(int * factor1, int factor1_l, int * factor2, int factor2_l, int * result, int * result_l){
    int num_zeros = 0, sum[factor1_l + factor2_l], sum_l = 0, tmp[factor1_l + factor2_l], tmp_l = 0;
    for(int i = factor2_l - 1; i >= 0; i--){
        multiply_function_val(factor1, factor1_l, factor2[i], tmp, &tmp_l);
        for(int j = 0; j < num_zeros; j++){
            tmp[tmp_l] = 0; 
            tmp_l++;
        }
        add_function(sum, sum_l, tmp, tmp_l, sum, &sum_l);
        num_zeros++;
    }
    *result_l = sum_l;
    for(int i = 0; i < sum_l; i++){
        result[i] = sum[i];
    }
}

void add_function(int * addend1, int addend1_l, int * addend2, int addend2_l, int * result, int * result_l){
    int max = addend1_l < addend2_l?addend2_l:addend1_l; 
    int r = 0, tmp_sum, i1 = addend1_l - 1, i2 = addend2_l - 1, rev_result[max + 1]; 
    int i = 0;
    while(i1 >= 0 && i2 >= 0){
        tmp_sum = addend1[i1] + addend2[i2] + r;
        rev_result[i] = tmp_sum % base;
        r = tmp_sum / base;
        i++;
        i1--;
        i2--;
    }
    if(i1 >= 0){
        while(i1 >= 0){
            tmp_sum = addend1[i1] + r;
            rev_result[i] = tmp_sum % base;
            r = tmp_sum / base;
            i1--;
            i++;
        }
    }
    if(i2 >= 0){
        while(i2 >= 0){
            tmp_sum = addend2[i2] + r;
            rev_result[i] = tmp_sum % base;
            r = tmp_sum / base;
            i2--;
            i++;
        }
    }
    if(r > 0){
        rev_result[i] = r;
        i++;
    }
    for(int k = 0; k < i; k++){
        result[k] = rev_result[i - 1 - k];
    }
    *result_l = i;
}