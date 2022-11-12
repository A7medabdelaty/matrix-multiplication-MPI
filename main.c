#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mpi.h"
#include <math.h>

#define N 4

// Print matrix function declaration
void printMatrix(int matrix[N][N]);

int matrix1[N][N];
int matrix2[N][N];
int productMatrix[N][N];


int main()
{
	int rows;
	int matrixSubset;

	int comm_sz;
	int my_rank;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	int numberOfWorkers = comm_sz - 1;

	if (my_rank == 0)
	{
		printf("\nMultiplying a %dx%d matrix using %d processor(s).\n\n", N, N, comm_sz);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				matrix1[i][j] = (rand() % 5) + 1;
				matrix2[i][j] = (rand() % 5) + 1;
			}
		}

		rows = ceil(N / numberOfWorkers);

		matrixSubset = 0;

		// Iterate through all of the workers and assign work
		for (int destinationProcessor = 1; destinationProcessor <= numberOfWorkers; destinationProcessor++)
		{
			MPI_Send(&matrixSubset, 1, MPI_INT, destinationProcessor, 1, MPI_COMM_WORLD);

			MPI_Send(&rows, 1, MPI_INT, destinationProcessor, 1, MPI_COMM_WORLD);

			MPI_Send(&matrix1[matrixSubset][0], rows * N, MPI_INT, destinationProcessor, 1, MPI_COMM_WORLD);

			MPI_Send(&matrix2, N * N, MPI_INT, destinationProcessor, 1, MPI_COMM_WORLD);

			matrixSubset = matrixSubset + rows;
		}

		// Retrieve results from all workers processors
		for (int sourceProcessor = 1; sourceProcessor <= numberOfWorkers; sourceProcessor++)
		{
			MPI_Recv(&matrixSubset, 1, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&rows, 1, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&productMatrix[matrixSubset][0], rows * N, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		printf("Matrix 1:\n");
		printMatrix(matrix1);
		printf("Matrix 2:\n");
		printMatrix(matrix2);
		printf("Product Matrix:\n");
		printMatrix(productMatrix);

	}
	else if (my_rank != 0)
	{
		int sourceProcessor = 0;
		MPI_Recv(&matrixSubset, 1, MPI_INT, sourceProcessor, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&rows, 1, MPI_INT, sourceProcessor, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&matrix1, rows * N, MPI_INT, sourceProcessor, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&matrix2, N * N, MPI_INT, sourceProcessor, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		/* Perform matrix multiplication */
		for (int k = 0; k < N; k++)
		{
			for (int i = 0; i < rows; i++)
			{
				productMatrix[i][k] = 0.0;
				for (int j = 0; j < N; j++)
				{
					productMatrix[i][k] = productMatrix[i][k] + matrix1[i][j] * matrix2[j][k];
				}
			}
		}

		MPI_Send(&matrixSubset, 1, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD);
		MPI_Send(&rows, 1, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD);
		MPI_Send(&productMatrix, rows * N, MPI_INT, sourceProcessor, 2, MPI_COMM_WORLD);
	}

	MPI_Finalize();
}

void printMatrix(int matrix[N][N])
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%d\t", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}