#include <iostream>

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "graphio.h"
#include "graph.h"
#include <stdbool.h>
}

#include <unordered_set>

char gfile[2048];

#define CHUNK_SIZE 500

void usage() {
	printf("./bfs <filename> <sourceIndex>\n");
	exit(0);
}


/*
You can ignore the ewgths and vwghts. They are there as the read function expects those values
row_ptr and col_ind are the CRS entities. nov is the Number of Vertices
*/

//void merge(struct node arr[], int l, int m, int r);
int FindNumberOfAdjacent(int current);
void LargestDegreeOrder(int nov);
//void mergeSort(struct node arr[], int l, int r);

etype *row_ptr;
vtype *col_ind;
int *color;
struct node *degree_array;
struct node *to_be_colored_array;

etype *d1_row_ptr;
vtype *d1_col_ind;

int* adjancency_count_array;

struct node
{
	int vertex_num;
	bool should_be_colored;

};

int FindMaxColor(int nov)
{
	int max = 0;

	for (int i = 0; i < nov; i++)
	{
		if (color[i] > max)
		{
			max = color[i];
		}
	}

	return max;

}

int FindNumberOfAdjacent(int current)
{
	int count = row_ptr[current + 1] - row_ptr[current];
	return count;

}

void LargestDegreeOrder(int nov)
{
	degree_array = new node[nov];

	int i;
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
	for (i = 0; i < nov; i++)
	{
		degree_array[i].vertex_num = i;
		degree_array[i].should_be_colored = true;
	}

}

int FindMaxAdj(int nov)
{
	int i;
	int largestdegree = 0;

#pragma omp parallel for reduction(max:largestdegree) schedule(dynamic, CHUNK_SIZE) 
	for (i = 0; i < nov; i++)
	{
		color[i] = 0;

		if (FindNumberOfAdjacent(i) > largestdegree)
		{
			largestdegree = FindNumberOfAdjacent(i);
		}

	}

	return largestdegree;
}

void ConvertToDistance1(int nov)
{
	std::unordered_set<int> * adjancency_map = new std::unordered_set<int>[nov];


	int i;
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
	for (i = 0; i < nov; i++)
	{
		for (int w = d1_row_ptr[i]; w < d1_row_ptr[i + 1]; w++)
		{
			int current_adj = d1_col_ind[w];

			//add distance-1 adj
			if (adjancency_map[i].find(current_adj) == adjancency_map[i].end())
			{
				adjancency_map[i].insert(current_adj);
			}

			for (int w2 = d1_row_ptr[current_adj]; w2 < d1_row_ptr[current_adj + 1]; w2++)
			{
				//add distance-2 adj
				int current_d2_adj = d1_col_ind[w2];

				if (adjancency_map[i].find(current_d2_adj) == adjancency_map[i].end())
				{
					adjancency_map[i].insert(current_d2_adj);
				}

			}

		}

		//delete self

		//for(std::unordered_set<int>::const_iterator it2 = adjancency_map[16].begin(); it2 != adjancency_map[16].end(); it2++)
		//{
		//	std::cout << *it2 << " ";
		//}

		//std::cout << std::endl;

		adjancency_map[i].erase(i);

	}

	int sum = 0;
	int j;

	row_ptr = new unsigned int[nov + 1];

	//#pragma omp parallel for reduction(+:sum) schedule(dynamic, CHUNK_SIZE)
	for (j = 0; j < nov; j++)
	{
		row_ptr[j] = sum;
		sum += adjancency_map[j].size();
	}
	row_ptr[nov] = sum;


	col_ind = new int[sum];


	std::unordered_set<int>::const_iterator it;

#pragma omp parallel for schedule(dynamic, CHUNK_SIZE) private(it)
	for (i = 0; i < nov; i++)
	{
		it = adjancency_map[i].begin();

		for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++)
		{
			col_ind[j] = *it;

			it++;
		}

	}

	delete[] adjancency_map;

}


int main(int argc, char *argv[])
{
	ewtype *ewghts;
	vwtype *vwghts;
	vtype nov, source;

	if (argc != 3)
		usage();

	const char* fname = argv[1];
	strcpy(gfile, fname);
	source = atoi(argv[2]);

	if (read_graph(gfile, &d1_row_ptr, &d1_col_ind, &ewghts, &vwghts, &nov, 0) == -1) {
		printf("error in graph read\n");
		exit(1);
	}

	for (int th = 0; th < 5; th++)
	{

		int NUM_THREADS = (2 << th) / 2;
		/****** YOUR CODE GOES HERE *******/
		omp_set_num_threads(NUM_THREADS);

		omp_set_dynamic(0);

		color = (int*)calloc(nov, sizeof(int));

		double start, end, prepfinish;

		start = omp_get_wtime();

		ConvertToDistance1(nov);

		LargestDegreeOrder(nov);

		prepfinish = omp_get_wtime();

		int largest_number_colors = FindMaxAdj(nov) + 1;

		int to_be_colored_number = nov;

		to_be_colored_array = degree_array;

		while (to_be_colored_number != 0)
		{
			int v;
			int current;
			int x;
			bool found;
			int * forbiddenColors;

#pragma omp parallel for schedule(dynamic, CHUNK_SIZE) private(current, x, found, forbiddenColors)
			for (v = 0; v < to_be_colored_number; v++)
			{
				current = to_be_colored_array[v].vertex_num;

				forbiddenColors = (int*)calloc(largest_number_colors, sizeof(int));


				for (int j = 1; j < largest_number_colors; j++)
				{

					forbiddenColors[j] = 0;
				}

				for (int w = row_ptr[current]; w < row_ptr[current + 1]; w++)
				{
					forbiddenColors[(color[col_ind[w]])] = 1;
				}

				found = false;

				x = 1;

				while (!found)
				{
					if (forbiddenColors[x] != 1)
					{
						found = true;
					}

					x++;
				}

				color[current] = x - 1;

				to_be_colored_array[v].should_be_colored = false;

				free(forbiddenColors);
			}

			///// to_be_colored_number - sum1 gives the number of remaining uncolored /////

			/////  conflict detection ///////

			int sum2 = 0;
#pragma omp parallel for reduction(+:sum2) schedule(dynamic, CHUNK_SIZE) private(current)
			for (int v = 0; v < to_be_colored_number; v++)
			{
				current = to_be_colored_array[v].vertex_num;

				for (int w = row_ptr[current]; w < row_ptr[current + 1]; w++)
				{
					if (color[current] == color[col_ind[w]] && current > col_ind[w])
					{
						to_be_colored_array[v].should_be_colored = true;

						sum2++;
						
						break;
					}
				}
			}
			//std::cout << "conflict number = " << sum2 << std::endl;

			int new_to_be_colored_number = sum2;

			//std::cout << "conflict number = " << new_to_be_colored_number << std::endl;

			///// U <- R /////////

			struct node * new_to_be_colored = (node*)calloc(new_to_be_colored_number, sizeof(struct node));
			int j = 0;
			for (v = 0; v < to_be_colored_number; v++)
			{
				if (to_be_colored_array[v].should_be_colored == true)
				{
					new_to_be_colored[j].vertex_num = to_be_colored_array[v].vertex_num;
					new_to_be_colored[j].should_be_colored = to_be_colored_array[v].should_be_colored;

					j++;
				}
			}

			free(to_be_colored_array);

			to_be_colored_number = new_to_be_colored_number;
			to_be_colored_array = new_to_be_colored;
      
      //std::cout << "to be colored number = " << to_be_colored_number << std::endl;


		}
		end = omp_get_wtime();

		//for (int i = 0; i < nov; i++)
		//{
		//	printf("%d ", degree_array[i].vertex_num);
		//}

		//printf("\n");

		//for (int i = 0; i < nov; i++)
		//{
		//	printf("%d ", degree_array[i].number_of_adj);
		//}

		//printf("\n");

		//for (int z = 0; z < nov; z++)
		//{
		// 		printf("%d ", color[z]);
		//}

		int sum3 = 0;

		for (int v = 0; v < nov; v++)
		{
			int current = v;

			for (int w = row_ptr[current]; w < row_ptr[current + 1]; w++)
			{
				if (color[current] == color[col_ind[w]])
				{
					sum3++;
				}
			}
		}

		printf("Number of Threads: %d ", (2 << th) / 2);

		printf("Final Conflict Amount: %d ", sum3);

		double elapsed_secs = end - start;

		printf("\nExecution Time: %f", elapsed_secs);

		printf("\nPreprocess Time: %f\n", prepfinish - start);

		printf("Max Color: %d\n\n\n", FindMaxColor(nov));


		free(color);
		degree_array = NULL;
		free(adjancency_count_array);
		free(row_ptr);
		free(col_ind);

	}

	free(d1_row_ptr);
	free(d1_col_ind);

	std::cin.get();


	return 1;
}
