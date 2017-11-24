# Parallelized Distance-2 Graph Coloring
This program executes sequential and parallel versions of distance-2 graph coloring algorithm, outputs the execution
times of sequential and parallel versions(with 1, 2, 4, 8, 16 threads) as well as the total number of color used in
each execution and runs a correctness check to confirm no distance-2 neighbors has the same color.

**Implemented by:**

 * @mbenlioglu
 * @ardaduz


## Getting Started
There are two Makefiles and a Python script included this project in order to crawl the data used for developement and
testing (Linux/Windows), and to compile sources into executable (Linux).

### Data Crawl

**For Linux:**

One of the Makefiles, `./Makefile-data`, is to crawl the data that is used for testing and developement. Running following
command on terminal will download and extract the data into `./data/`. Note that, this folder will be created if it does not
already exists.

    $ make -f Makefile-data

Also note that, after crawling and extracting data, following command can be run to get rid of downloaded archive files (Not
the extracted data folders/files just archives).

    $ make -f Makefile-data clean

**For Windows:**

Since Makefiles and command-line tools used in the Makeflie (e.g. tar) is not officially supported in Windows, we have
provided a script written in [Python 2.7.x](https://docs.python.org/2/), to download and extract the data. Note that,
the data is archived with TAR & GZ compression, where there is also no official tools from Windows to decompress.
Therefore, we have used [WinRAR](https://www.win-rar.com/) to decompress these archieves in the script, which you may
need to download if you do not have. Running the following script will download and extract the data into `./data/` folder,
and will create this folder if it doesn't already exists.
    
    $ python ./DataCrawl-windows.py

Note that, running this script with no parameters will assume that, `"C:/Program Files/WinRAR/WinRAR.exe"` exists. If
WinRAR is installed in some other location, it can be specified with `--winrar_path` flag. Also note that giving `-r` or
`--remove-archieves` flag will remove original archives after donwloading and extraction finishes (similar to `make clean`)_.
For more information about this script, execute following command in project root.

    $ python ./DataCrawl-windows.py -h
------
Input data is in [MTX format](http://math.nist.gov/MatrixMarket/formats.html) and converted into [CRS format](http://netlib.org/linalg/html_templates/node91.html)
for processing. For more information about the data used in this project visit [SuiteSparse Matrix Collection](https://sparse.tamu.edu/).

### Compiling

The other Makefile, `./Makefile` is to compile codes into executable in Linux environment. Executing following command
in the project root will create an executable named "coloring" which takes single input argument, which is the path to
the dataset.

    $ make

For Windows use Visual Studio with provided *.sln & *.vxproj files and use /openmp flag to compile the codes.

The program will run sequential and parallel versions of the algorithm and output the execution time and number of colors used.


## Algorithm & Implementation
In this section, sequential and parallel version of the algorithm will be explained. Performance results of each algorithm can
be found in results section.

Note that, colors are represented as short integers starting from 0 to 32,767. This is because, we know that number of maximum
colors used in a graph is bounded by number of maximum distance-2 neighbours a vertex can have, and since we are using sparse
graphs, this number doesn’t even get close to the limit.

### Direct Approach
#### Sequential Version
Since the size of the graphs are huge, the algorithm follows a greedy approach instead of dynamic programming, therefore total
number of colors used in the graph may result in a suboptimal number. In this approach, all the vertices are iterated over one-by-one
to find appropriate color that satisfies distance-1 coloring conditions. 

On each iteration (i.e. for every vertex), the smallest color that is not used by its neighboring vertices is assigned to the vertex.
Therefore, this implantation favors a first-fit manner for assigned colors.

In order to find the smallest available color, neighbors of the vertex are visited, and a Boolean array is used to store this
information. Then, by iterating over this Boolean smallest unused color can be easily determined. By this approach, we can
determine smallest unused color in O(k) time complexity, where k represents the number of neighbors of the vertex.

#### Parallel version
Parallelized version, uses the same approach as the sequential version, only difference is that multiple vertices are controlled
and assigned at the same time, which introduces a race condition. Therefore, at the end of initial colorization step, a conflict
detection phase occurs, where all graph is checked for conflicted vertices (i.e. distance-2 neighbors with same color) and these
vertices recolored in parallel. Conflict detection phase repeats until there are no conflicts left on the graph.

## Results
In this section, result tables containing various information about execution results is given.

### Direct Approach

![Execution Time Result Table](/res/time_results_direct.png)

**Table 1.1.:** Executions times, speedups and efficiencies of implementations explained in the "Direct Approach" section.

![Color/Conflict fix Table](/res/color_conf_results_direct.png)

**Table 1.2.:** Number of colors used, and number of conflict fix iterations of implementations explained in the "Direct Approach"
section.

