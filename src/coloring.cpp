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
#include<unordered_set>
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
void distance2ToDistance1(etype *rowD1, vtype *colD1, vtype nov, etype *rowD2, vtype *colD2)
{
	std::unordered_set<int> * adjancency_map = new std::unordered_set<int>[nov];


	int i;
#pragma omp parallel for
	for (i = 0; i < nov; i++)
	{
		for (int w = rowD1[i]; w < rowD1[i + 1]; w++)
		{
			int current_adj = colD1[w];

			//add distance-1 adj
			if (adjancency_map[i].find(current_adj) == adjancency_map[i].end())
			{
				adjancency_map[i].insert(current_adj);
			}

			for (int w2 = rowD1[current_adj]; w2 < rowD1[current_adj + 1]; w2++)
			{
				//add distance-2 adj
				int current_d2_adj = colD1[w2];

				if (adjancency_map[i].find(current_d2_adj) == adjancency_map[i].end())
				{
					adjancency_map[i].insert(current_d2_adj);
				}

			}

		}

		//delete self

		adjancency_map[i].erase(i);

	}

	int sum = 0;
	int j;

	rowD2 = new unsigned int[nov + 1];

	for (j = 0; j < nov; j++)
	{
		rowD2[j] = sum;
		sum += adjancency_map[j].size();
	}
	rowD2[nov] = sum;


	colD2 = new int[sum];


	std::unordered_set<int>::const_iterator it;

#pragma omp parallel for private(it)
	for (i = 0; i < nov; i++)
	{
		it = adjancency_map[i].begin();

		for (int j = rowD2[i]; j < rowD2[i + 1]; j++)
		{
			colD2[j] = *it;

			it++;
		}

	}

	delete[] adjancency_map;

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
		int c, colStart, colEnd, conflictIndex, temp;
		int i, j, k;
#pragma omp parallel for private(j, k, c, colStart, colEnd, conflictIndex, temp)
		for (i = 0; i < nov; i++)
		{
			c = colors[i];
			colStart = row[i];
			colEnd = row[i + 1];
			for (j = colStart; j < colEnd; j++)
			{
				if (colors[col[j]] == c)
				{
					conflictIndex = i < col[j] ? i : col[j];
					if (!isDetected[conflictIndex])
					{
						isDetected[conflictIndex] = true;
#pragma omp atomic capture
						temp = index++;
						out[temp] = conflictIndex;
					}
				}
			}
		}

		// reset isDetected array
#pragma omp parallel for 
		for (int i = 0; i < index; i++)
			isDetected[out[i]] = false;

		return index;
	}

	namespace Direct
	{
		void resetIsColorUsed(int vertex, etype *row, vtype *col, int colors[], bool isColorUsed[])
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
		}
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
				{
					resetIsColorUsed(vertex, row, col, colors, isColorUsed);
					return i;
				}

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
			int i, c;
			startTime = omp_get_wtime();
#pragma omp parallel for private(c)
			for (i = 0; i < nov; i++)
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
#pragma omp for private(c)
				for (i = 0; i < confCnt; i++) // recolor
				{
					c = getSmallestAvailableColor(conflictedVertices[i], row, col, nov, colors, isColorUsed, maxEdgeCnt);
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
		int c, colStart, colEnd, d2colStart, d2colEnd, conflictIndex, temp;
		int i, j, k;
#pragma omp parallel for private(j, k, c, colStart, colEnd, d2colStart, d2colEnd, conflictIndex, temp)
		for (i = 0; i < nov; i++)
		{
			c = colors[i];
			colStart = row[i];
			colEnd = row[i + 1];
			for (j = colStart; j < colEnd; j++)
			{
				if (colors[col[j]] == c)
				{
					conflictIndex = i < col[j] ? i : col[j];
					if (!isDetected[conflictIndex])
					{
						isDetected[conflictIndex] = true;
#pragma omp atomic capture
						temp = index++;
						out[temp] = conflictIndex;
					}
				}

				d2colStart = row[col[j]];
				d2colEnd = row[col[j] + 1];
				for (k = d2colStart; k < d2colEnd; ++k)
				{
					if (colors[col[k]] == c && col[k] != i)
					{
						conflictIndex = i < col[k] ? i : col[k];
						if (!isDetected[conflictIndex])
						{
							isDetected[conflictIndex] = true;
#pragma omp atomic capture
							temp = index++;
							out[temp] = conflictIndex;
						}
					}
				}
			}
		}

		// reset isDetected array
#pragma omp parallel for 
		for (i = 0; i < index; i++)
			isDetected[out[i]] = false;

		return index;
	}

	namespace Direct
	{
		void resetIsColorUsed(int vertex, etype *row, vtype *col, int colors[], bool isColorUsed[])
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
		}

		int getSmallestAvailableColor(int vertex, etype *row, vtype *col, vtype nov, int colors[], bool isColorUsed[])
		{
			int colStart, colEnd;
			colStart = row[vertex];
			colEnd = row[vertex + 1];

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
				{
					resetIsColorUsed(vertex, row, col, colors, isColorUsed);
					return i;
				}

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
				isColorUsed = new bool[nov + 1]();
			}

			// first stage coloring
			int i, c;
			startTime = omp_get_wtime();
#pragma omp parallel for private(c)
			for (i = 0; i < nov; i++)
			{
				c = getSmallestAvailableColor(i, row, col, nov, colors, isColorUsed);
				colors[i] = c;
			}

			int confCnt = 0;
			// fix conflict
			do
			{
				// detect conflicted vertices
				confCnt = detect_conflicts(row, col, nov, colors, isVertexDetected, conflictedVertices);
#pragma omp for private(c)
				for (i = 0; i < confCnt; i++) // recolor
				{
					c = getSmallestAvailableColor(conflictedVertices[i], row, col, nov, colors, isColorUsed);
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
}
//===========================================================================================================================
//===========================================================================================================================

int main(int argc, char *argv[])
{
	std::string baseName = std::string(argv[0]);
	std::string fillerAsterisk(100, '*');
	std::string fillerDashes(100, '-');

	// Get executable name from path
#ifdef _WIN32
	baseName = baseName.substr(baseName.rfind('\\') + 1);
#else
	baseName = baseName.substr(baseName.rfind('/') + 1);
#endif // _WIN32

	std::cout << fillerAsterisk << std::endl;
	std::cout << "Starting " << baseName << "...\n";
	std::cout << fillerAsterisk << std::endl;

	etype *row_ptr;
	vtype *col_ind;
	ewtype *ewghts;
	vwtype *vwghts;
	vtype nov;

	std::cout << "\nReading graph... ";
	// Graph reading
	if (argc < 2)
	{
		print_usage();
		return EXIT_FAILURE;
	}
	
	if (read_graph(argv[1], &row_ptr, &col_ind, &ewghts, &vwghts, &nov, 0) == -1)
	{
		std::cout << "error in graph read\n";
		return EXIT_FAILURE;
	}
	std::cout << "done\n\n";

	// Analyse graph to find maximum edge count on a vertex
	int maxEdgeCnt = findMaxEdgeCnt(row_ptr, col_ind, nov);

	// Performance analysis
	perfData perfSeq, perfPar[5];

	int *colors = new int[nov];
	std::fill_n(colors, nov, -1);

	//===========================================================================================================================
	// Direct Approach
	//===========================================================================================================================
	std::cout << fillerDashes << "\n";
	std::cout << "Starting performance analysis for direct approch\n";
	std::cout << fillerDashes << "\n";
	
	// Sequential
	std::cout << "Starting sequential algorithm...";
	std::cout.flush();
	perfSeq = Distance2::Direct::color_graph_seq(row_ptr, col_ind, nov, colors);
	std::cout << "ended\n";
#ifdef DEBUG
	int outSize = nov;

	// these two are used in the detect_conflicts, for correctness we only need to check conflict count.
	bool *isDetected = new bool[outSize]();
	int *out = new int[outSize]();

	std::cout << "Running correctness check...";
	std::cout.flush();
	omp_set_num_threads(1);
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
		std::cout.flush();
		perfPar[i] = Distance2::Direct::color_graph_par(row_ptr, col_ind, nov, colors);
		std::cout << "ended\n";
#ifdef DEBUG
		std::cout << "Running correctness check...";
		std::cout.flush();
		omp_set_num_threads(1);
		s = !Distance2::detect_conflicts(row_ptr, col_ind, nov, colors, isDetected, out) ? "correct\n" : "wrong!\n";
		std::cout << s;
		std::fill_n(isDetected, outSize, false);
#endif // DEBUG
		std::fill_n(colors, nov, -1); // reinitialize
	}

	// Print results
	std::cout << fillerDashes << "\n";
	printf("| %-15s | %-12s | %-15s | %-13s | %-16s |\n", "Algorithm", "Thread Count", "# of Conf.Fixes", "# of Colors", "Exec. Time");
	std::cout << fillerDashes << "\n";
	printf("| %-15s | %-12d | %-15d | %-13d | %-14.10f s |\n",
		"Sequential", 1, perfSeq.mergeConflictCnt, perfSeq.colorCnt, perfSeq.execTime);
	for (int i = 0; i < 5; i++)
		printf("| %-15s | %-12d | %-15d | %-13d | %-14.10f s |\n",
			"Parallel", (2 << i) / 2, perfPar[i].mergeConflictCnt, perfPar[i].colorCnt, perfPar[i].execTime);
	std::cout << "\n";

	//===========================================================================================================================
	// D2-D1 conversion then D1 solution
	//===========================================================================================================================
	
	std::cout << fillerDashes << "\n";
	std::cout << "Starting performance analysis for D2-to-D1 conversion approch\n";
	std::cout << fillerDashes << "\n";

	etype *d2_row_ptr = NULL;
	vtype *d2_col_ind = NULL;
	double startTime, prepTime;
	double parPrepTime[5];

	// Preprocessing: convert graph to d1 graph
	std::cout << "Preprocessing...";
	std::cout.flush();
	startTime = omp_get_wtime();
	distance2ToDistance1(row_ptr, col_ind, nov, d2_row_ptr, d2_col_ind);
	prepTime = omp_get_wtime() - startTime;
	std::cout << " done\n";

	// Sequential
	std::cout << "Starting sequential algorithm...";
	std::cout.flush();
	perfSeq = Distance1::Direct::color_graph_seq(d2_row_ptr, d2_col_ind, nov, colors, maxEdgeCnt);
	std::cout << "ended\n";
#ifdef DEBUG
	// these two are used in the detect_conflicts, for correctness we only need to check conflict count.
	isDetected = new bool[outSize]();
	out = new int[outSize]();

	std::cout << "Running correctness check...";
	std::cout.flush();
	s = !Distance1::detect_conflicts(d2_row_ptr, d2_col_ind, nov, colors, isDetected, out) ? "correct\n" : "wrong!\n";
	std::cout << s;
	std::fill_n(isDetected, outSize, false);
#endif // DEBUG
	std::fill_n(colors, nov, -1); // reinitialize
	delete[] d2_row_ptr;
	delete[] d2_col_ind;

	// Parallel
	for (int i = 0; i < 5; i++)
	{
		omp_set_num_threads((2 << i) / 2);

		// Preprocessing: convert graph to d1 graph
		std::cout << "Preprocessing...";
		std::cout.flush();
		startTime = omp_get_wtime();
		distance2ToDistance1(row_ptr, col_ind, nov, d2_row_ptr, d2_col_ind);
		parPrepTime[i - 1] = omp_get_wtime() - startTime;
		std::cout << " done\n";

		std::cout << "Starting parallel algorithm with " << (2 << i) / 2 << " threads...";
		std::cout.flush();
		perfPar[i - 1] = Distance1::Direct::color_graph_par(d2_row_ptr, d2_col_ind, nov, colors, maxEdgeCnt);
		std::cout << "ended\n";
#ifdef DEBUG
		std::cout << "Running correctness check...";
		std::cout.flush();
		s = !Distance1::detect_conflicts(d2_row_ptr, d2_col_ind, nov, colors, isDetected, out) ? "correct\n" : "wrong!\n";
		std::cout << s;
		std::fill_n(isDetected, outSize, false);
#endif // DEBUG
		std::fill_n(colors, nov, -1); // reinitialize
	}

	// Print results
	printf("| %-15s | %-12s | %-15s | %-12s | %-15s | %-15s |\n", "Algorithm", "Thread Count", "# of Conf.Fixes", "# of Colors", "Prep. Time", "Exec. Time");
	std::cout << fillerDashes << "\n";
	printf("| %-15s | %-12d | %-15d | %-12d | %-12.10f s | %-12.10f s |\n",
		"Sequential", 1, perfSeq.mergeConflictCnt, perfSeq.colorCnt, prepTime, perfSeq.execTime);
	for (int i = 1; i < 5; i++)
		printf("| %-15s | %-12d | %-15d | %-12d | %-12.10f s | %-12.10f s |\n",
			"Parallel", (2 << i) / 2, perfPar[i].mergeConflictCnt, perfPar[i].colorCnt, parPrepTime[i], perfPar[i].execTime);
	std::cout << "\n";



	std::cout << "Press Enter to continue...\n";
	std::cin.get();

	return 0;
}