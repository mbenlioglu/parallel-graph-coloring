/*
* Created by: mbenlioglu, October 26, 2017
*
* This program executes sequential and parallel versions of distance-1 graph coloring algorithm, outputs the execution
* times of sequential and parallel versions(with 1, 2, 4, 8, 16 threads) as well as the total number of color used in
* each execution and runs a correctness check to confirm no adjacent vertices has the same color.
*
*/

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graphio.h"
#include "graph.h"
#ifdef __cplusplus
}
#endif /*__cplusplus*/

#include <iostream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <omp.h>


#define DEBUG

//===========================================================================================================================
// Common Structs/Enums
//===========================================================================================================================
/*
	Struct for reporting performance analysis results
*/
typedef	struct perfData
{
	double prepTime, execTime;
	int colorCnt, mergeConflictCnt;
}perfData;

//===========================================================================================================================
// Common Functions
//===========================================================================================================================


/*
	Finds and returns total number of colors used in the coloring
*/
inline int num_of_colors(vtype nov, int colors[])
{
	int largest = -1;

	for (int i = 0; i < nov; i++)
		if (colors[i] > largest)
			largest = colors[i];

	return largest + 1;
}

/*
	Finds and returns the maximum number of neighbours a vertex can have in the given graph
*/
int findMaxEdgeCnt(etype *row, vtype *col, vtype nov)
{
	int max = -1;

	for (int i = 0; i < nov; ++i)
	{
		int edgeCnt = row[i + 1] - row[i];
		if (edgeCnt > max)
			max = edgeCnt;
	}

	return max;
}

/*
	Update *row, *col pointers to make distance-2 nighbors, distance-1 neighbors
*/
void distance2ToDistance1(etype *row, vtype *col, vtype nov)
{
	
}

void print_usage()
{
	std::cout << "./coloring <filename>\n";
	std::cout << "Press Enter to continue...";
	std::cin.get();
}

//===========================================================================================================================
// Namespaces for implementations
//===========================================================================================================================

// Distance-1 graph coloring
namespace Distance1
{
	/*
	Traverses entire graph to find conflicts (i.e. adjecent vertices with same color), returns the number of such
	conflicts as return value, indices of these vertices are written into out array.

	*row and *col: pointers define the starting point of the graph
	nov:		   number of vertices in the graph
	colors:		   the array storing color values assigned to each vertex (-1 means unassigned)
	isDetected:    the array that will be used as temporal storage to find whether the vertex is already marked as
	conflict, initially all values assumed to be false (taken as parameter because of efficiency considerations)
	out:		   output array that contains all the vertices marked to be re-colored (size >= nov/2)

	returns:       number of conflicts detected.
	*/
	int detect_conflicts(etype *row, vtype *col, vtype nov, int colors[], bool isDetected[], int out[])
	{
		unsigned int index = 0;

		for (int i = 0; i < nov; i++)
		{
			int c = colors[i];
			int colStart = row[i], colEnd = row[i + 1];
			for (int j = colStart; j < colEnd; j++)
			{
				if (colors[col[j]] == c)
				{
					int conflictIndex = i < col[j] ? i : col[j];
					if (!isDetected[conflictIndex])
					{
						isDetected[conflictIndex] = true;
						out[index++] = conflictIndex;
					}
				}
			}
		}

		// reset isDetected array
		for (int i = 0; i < index; i++)
			isDetected[out[i]] = false;

		return index;
	}

	namespace Direct
	{
		/*
			Traverses the neighbours of the given vertex and returns the smallest available color for the vertex

			*row and *col: pointers define the starting point of the graph
			vertex:		   the index of the vertex that wanted to be colored
			nov:		   number of vertices in the graph
			colors:		   the array storing color values assigned to each vertex (-1 means unassigned)
			isColorUsed:   the array that will be used as temporal storage to find which is the smallest unused color value
				(taken as parameter because of efficiency considerations)
			maxEdgeCnt:    maximum number of edges a vertex can have (i.e. size of isColorUsed - 1)
		*/
		int getSmallestAvailableColor(int vertex, etype *row, vtype *col, vtype nov, int colors[], bool isColorUsed[], int maxEdgeCnt)
		{
			int colStart, colEnd;
			colStart = row[vertex];
			colEnd = row[vertex + 1];

			//reset setted parts in bool array
			for (int j = colStart; j < colEnd; j++)
			{
				int c = colors[col[j]];
				if (c >= 0)
					isColorUsed[c] = false;
			}

			// track whether a color is used it not
			for (int i = colStart; i < colEnd; i++)
			{
				int c = colors[col[i]];
				if (c >= 0)
					isColorUsed[c] = true;
			}

			//	return the smallest unused color
			for (int i = 0; i < maxEdgeCnt + 1; i++)
				if (!isColorUsed[i])
					return i;

			// should never enter here
			throw std::runtime_error("Should NEVER enter HERE!!");
		}

		perfData color_graph_seq(etype *row, vtype *col, vtype nov, int colors[], int maxEdgeCnt)
		{
			perfData result;
			double startTime, endTime;
			int maxColor = 0;
			bool *isColorUsed = new bool[maxEdgeCnt + 1]();

			startTime = omp_get_wtime();
			for (int i = 0; i < nov; i++)
			{
				int c = getSmallestAvailableColor(i, row, col, nov, colors, isColorUsed, maxEdgeCnt);
				colors[i] = c;
				if (c > maxColor) maxColor = c;
			}
			endTime = omp_get_wtime();

			// clean up
			delete[] isColorUsed;

			result.colorCnt = maxColor + 1;
			result.execTime = endTime - startTime;
			result.prepTime = result.mergeConflictCnt = 0;
			return result;
		}

		perfData color_graph_par(etype *row, vtype *col, vtype nov, int colors[], int maxEdgeCnt)
		{
			perfData result;
			double startTime, endTime;

			int mergeConflictCnt = -1;

			int confArrSize = nov / 2 + 1;
			int *conflictedVertices = new int[confArrSize]();
			bool *isVertexDetected = new bool[nov]();
			static bool *isColorUsed;
#pragma omp threadprivate(isColorUsed)

#pragma omp parallel
			{
				isColorUsed = new bool[maxEdgeCnt + 1]();
			}

			// first stage coloring
			startTime = omp_get_wtime();
#pragma omp parallel for
			for (int i = 0; i < nov; i++)
			{
				int c = getSmallestAvailableColor(i, row, col, nov, colors, isColorUsed, maxEdgeCnt);
				colors[i] = c;
			}

			int confCnt = 0;
			// fix conflict
			do
			{
				// detect conflicted vertices
				confCnt = detect_conflicts(row, col, nov, colors, isVertexDetected, conflictedVertices);
#pragma omp for
				for (int i = 0; i < confCnt; i++) // recolor
				{
					int c = getSmallestAvailableColor(conflictedVertices[i], row, col, nov, colors, isColorUsed, maxEdgeCnt);
					colors[conflictedVertices[i]] = c;
				}
				++mergeConflictCnt;
			} while (confCnt > 0);
			endTime = omp_get_wtime();

			// clean up
			delete[] isVertexDetected;
			delete[] conflictedVertices;
#pragma omp parallel
			{
				delete[] isColorUsed;
			}

			result.colorCnt = num_of_colors(nov, colors);
			result.execTime = endTime - startTime;
			result.mergeConflictCnt = mergeConflictCnt;
			result.prepTime = 0;
			return result;
		}
	}

	namespace Heuristic
	{
		perfData color_graph_seq(etype *row, vtype *col, vtype nov, int colors[])
		{
			throw std::runtime_error("Not Implemented");
		}

		perfData color_graph_par(etype *row, vtype *col, vtype nov, int colors[])
		{
			throw std::runtime_error("Not Implemented");
		}
	}
}

// Distance-2 Graph coloring
namespace Distance2
{
	int detect_conflicts(etype *row, vtype *col, vtype nov, int colors[], bool isDetected[], int out[])
	{
		unsigned int index = 0;

		for (int i = 0; i < nov; i++)
		{
			int c = colors[i];
			int colStart = row[i], colEnd = row[i + 1];
			for (int j = colStart; j < colEnd; j++)
			{
				if (colors[col[j]] == c)
				{
					int conflictIndex = i < col[j] ? i : col[j];
					if (!isDetected[conflictIndex])
					{
						isDetected[conflictIndex] = true;
						out[index++] = conflictIndex;
					}
				}

				int d2colStart = row[col[j]], d2colEnd = row[col[j] + 1];
				for (int k = d2colStart; k < d2colEnd; ++k)
				{
					if (colors[col[k]] == c && col[k] != i)
					{
						int conflictIndex = i < col[k] ? i : col[k];
						if (!isDetected[conflictIndex])
						{
							isDetected[conflictIndex] = true;
							out[index++] = conflictIndex;
						}
					}
				}
			}
		}

		// reset isDetected array
		for (int i = 0; i < index; i++)
			isDetected[out[i]] = false;

		return index;
	}

	namespace Direct
	{
		int getSmallestAvailableColor(int vertex, etype *row, vtype *col, vtype nov, int colors[], bool isColorUsed[])
		{
			int colStart, colEnd;
			colStart = row[vertex];
			colEnd = row[vertex + 1];

			//reset setted parts in bool array
			for (int i = colStart; i < colEnd; i++)
			{
				int c = colors[col[i]];
				if (c >= 0)
					isColorUsed[c] = false;

				int d2colStart = row[col[i]], d2colEnd = row[col[i] + 1];
				for (int j = d2colStart; j < d2colEnd; j++)
				{
					c = colors[col[j]];
					if (c >= 0 && col[j] != vertex)
						isColorUsed[c] = false;
				}
			}

			// track whether a color is used it not
			for (int i = colStart; i < colEnd; i++)
			{
				int c = colors[col[i]];
				if (c >= 0)
					isColorUsed[c] = true;

				int d2colStart = row[col[i]], d2colEnd = row[col[i] + 1];
				for (int j = d2colStart; j < d2colEnd; j++)
				{
					c = colors[col[j]];
					if (c >= 0 && col[j] != vertex)
						isColorUsed[c] = true;
				}
			}

			//	return the smallest unused color
			for (int i = 0; i < nov + 1; i++)
				if (!isColorUsed[i])
					return i;

			// should never enter here
			throw std::runtime_error("Should NEVER enter HERE!!");
		}

		perfData color_graph_seq(etype *row, vtype *col, vtype nov, int colors[])
		{
			perfData result;
			double startTime, endTime;
			int maxColor = 0;
			bool *isColorUsed = new bool[nov + 1]();

			startTime = omp_get_wtime();
			for (int i = 0; i < nov; i++)
			{
				int c = getSmallestAvailableColor(i, row, col, nov, colors, isColorUsed);
				colors[i] = c;
				if (c > maxColor) maxColor = c;
			}
			endTime = omp_get_wtime();

			// clean up
			delete[] isColorUsed;

			result.colorCnt = maxColor + 1;
			result.execTime = endTime - startTime;
			result.prepTime = result.mergeConflictCnt = 0;
			return result;
		}

		perfData color_graph_par(etype *row, vtype *col, vtype nov, int colors[])
		{
			// TODO:
			return perfData();
		}
	}
}
//===========================================================================================================================
//===========================================================================================================================

int main(int argc, char *argv[])
{
	etype *row_ptr;
	vtype *col_ind;
	ewtype *ewghts;
	vwtype *vwghts;
	vtype nov;

	// Graph reading
	if (argc < 2)
	{
		print_usage();
		return 1;
	}
	
	if (read_graph(argv[1], &row_ptr, &col_ind, &ewghts, &vwghts, &nov, 0) == -1)
	{
		std::cout << "error in graph read\n";
		return 1;
	}

	// Analyse graph to find maximum edge count on a vertex
	int maxEdgeCnt = findMaxEdgeCnt(row_ptr, col_ind, nov);

	// Performance analysis
	perfData perfSeq, perfPar[5];

	int *colors = new int[nov];
	std::fill_n(colors, nov, -1);

	//===========================================================================================================================
	// Direct Approach
	//===========================================================================================================================
	std::cout << std::setfill('*') << std::setw(100) << "-\n";
	std::cout << "Starting performance analysis for direct approch\n\n";
	std::cout << std::setfill('*') << std::setw(100) << "-\n";
	
	// Sequential
	std::cout << "Starting sequential algorithm...";
	perfSeq = Distance2::Direct::color_graph_seq(row_ptr, col_ind, nov, colors);
	std::cout << "ended\n";
#ifdef DEBUG
	int outSize = nov;

	// these two are used in the detect_conflicts, for correctness we only need to check conflict count.
	bool *isDetected = new bool[outSize]();
	int *out = new int[outSize]();

	std::cout << "Running correctness check...";
	std::string s = !Distance2::detect_conflicts(row_ptr, col_ind, nov, colors, isDetected, out) ? "correct\n" : "wrong!\n";
	std::cout << s;
	std::fill_n(isDetected, outSize, false);
#endif // DEBUG
	std::fill_n(colors, nov, -1); // reinitialize
	
	// Parallel
	for (int i = 0; i < 5; i++)
	{
		omp_set_num_threads((2 << i) / 2);
		std::cout << "Starting parallel algorithm with " << (2 << i) / 2 << " threads...";
		perfPar[i] = Distance2::Direct::color_graph_par(row_ptr, col_ind, nov, colors);
		std::cout << "ended\n";
#ifdef DEBUG
		std::cout << "Running correctness check...";
		s = !Distance2::detect_conflicts(row_ptr, col_ind, nov, colors, isDetected, out) ? "correct\n" : "wrong!\n";
		std::cout << s;
		std::fill_n(isDetected, outSize, false);
#endif // DEBUG
		std::fill_n(colors, nov, -1); // reinitialize
	}

	// Print results
	printf("| %-15s | %-12s | %-15s | %-12s | %-15s |\n", "Algorithm", "Thread Count", "# of Conf.Fixes", "# of Colors", "Exec. Time");
	std::cout << std::setfill('-') << std::setw(85) << "-\n";
	printf("| %-15s | %-12d | %-15d | %-12d | %-12.10f s |\n",
		"Sequential", 1, perfSeq.mergeConflictCnt, perfSeq.colorCnt, perfSeq.execTime);
	for (int i = 0; i < 5; i++)
		printf("| %-15s | %-12d | %-15d | %-12d | %-12.10f s |\n",
			"Parallel", (2 << i) / 2, perfPar[i].mergeConflictCnt, perfPar[i].colorCnt, perfPar[i].execTime);
	std::cout << "\n";

	//===========================================================================================================================
	// D2-D1 conversion then D1 solution
	//===========================================================================================================================
	





	std::cout << "Press Enter to continue...\n";
	std::cin.get();

	return 0;
}