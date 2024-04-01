#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define MAX_THREADS 100 // Maximum number of threads

// Structure to hold the data required by each thread
typedef struct {
    int row, col;       // Row and column index of the cell to compute
    int m, n, k;        // Dimensions of matrices
    double **A, **B, **C;  // Matrices A, B, and C
} ThreadData;

// Function prototypes
void *compute_cell(void *thread_arg);
void read_matrices_from_file(char *filename, double ***A, double ***B, int *m, int *k, int *n);
void free_matrix(double **matrix, int rows);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *filename = argv[1];

    int m, n, k; // Dimensions of matrices
    double **A, **B, **C;

    // Read matrices from the file
    read_matrices_from_file(filename, &A, &B, &m, &k, &n);

    // Allocate memory for matrix C
    C = (double **)malloc(m * sizeof(double *));
    for (int i = 0; i < m; i++)
        C[i] = (double *)calloc(n, sizeof(double));

    // Initialize pthread variables
    pthread_t threads[MAX_THREADS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Create threads for each cell of matrix multiplication
    int thread_count = 0;
    ThreadData thread_data[m][n];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            thread_data[i][j].row = i;
            thread_data[i][j].col = j;
            thread_data[i][j].m = m;
            thread_data[i][j].n = n;
            thread_data[i][j].k = k;
            thread_data[i][j].A = A;
            thread_data[i][j].B = B;
            thread_data[i][j].C = C;
            pthread_create(&threads[thread_count], &attr, compute_cell, (void *)&thread_data[i][j]);
            thread_count++;
        }
    }

    // Join threads
    for (int i = 0; i < m * n; i++)
        pthread_join(threads[i], NULL);

    // Display the result matrix C
    printf("Matrix C (Product of A and B):\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            printf("%.6lf ", C[i][j]);
        printf("\n");
    }

    // Free allocated memory
    for (int i = 0; i < m; i++) {
        free(A[i]);
        free(C[i]);
    }
    for (int i = 0; i < k; i++)
        free(B[i]);
    free(A);
    free(B);
    free(C);

    pthread_attr_destroy(&attr);
    pthread_exit(NULL);

    return 0;
}

// Function to compute a cell of matrix multiplication
void *compute_cell(void *thread_arg) {
    ThreadData *data = (ThreadData *)thread_arg;

    // Extract data from the argument
    int row = data->row;
    int col = data->col;
    int k = data->k;
    double **A = data->A;
    double **B = data->B;
    double **C = data->C;

    // Compute the value of the cell (row, col) in matrix C
    for (int i = 0; i < k; i++) {
        C[row][col] += A[row][i] * B[i][col];
    }

    pthread_exit(NULL);
}

// Function to read matrices from the file
void read_matrices_from_file(char *filename, double ***A, double ***B, int *m, int *k, int *n) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Read dimensions
    if (fscanf(file, "%d %d %d", m, k, n) != 3 || *m <= 0 || *k <= 0 || *n <= 0) {
        fprintf(stderr, "Invalid dimensions in the input file.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for matrices A and B
    *A = (double **)malloc((*m) * sizeof(double *));
    *B = (double **)malloc((*k) * sizeof(double *));
    if (*A == NULL || *B == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < *m; i++)
        (*A)[i] = (double *)malloc((*k) * sizeof(double));
    for (int i = 0; i < *k; i++)
        (*B)[i] = (double *)malloc((*n) * sizeof(double));
    if ((*A)[*m - 1] == NULL || (*B)[*k - 1] == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Read matrix A
    for (int i = 0; i < *m; i++) {
        for (int j = 0; j < *k; j++) {
            if (fscanf(file, "%lf", &((*A)[i][j])) != 1) {
                fprintf(stderr, "Error reading matrix A from the input file.\n");
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Read matrix B
    for (int i = 0; i < *k; i++) {
        for (int j = 0; j < *n; j++) {
            if (fscanf(file, "%lf", &((*B)[i][j])) != 1) {
                fprintf(stderr, "Error reading matrix B from the input file.\n");
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }
    }

    fclose(file);
}

// Function to free memory allocated for a matrix
void free_matrix(double **matrix, int rows) {
    for (int i = 0; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}
