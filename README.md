# IMPC: Influence Maximization based on Multi-Neighbor Potential in Community Networks
The Influence Maximization Algorithm based on Multi-neighbor Potential
Welcome to use the source code of the CoFIM algorithm proposed by Shang et al.

When using this code please cite our paper titled "IMPC: Influence Maximization based on Multi-Neighbor Potential in Community Networks" (to appear in the future).


- Compile

	Makefile is included, just type type "make" to compile all the source codes.
	"gcc 4.7.2" preferred


- Execute

	Example:
		./IMPC -data NetHEPT.txt -com NetHEPT_com.txt -k 50 -beta 0.2 -alpha 0.2

	Arguments:
		-data:
			the graph file
		-com:
			the community file
		-k:
			number of seed nodes
		-beta:
			the algorithm parameter
		-alpha:
			the algorithm parameter


- Evalution

	To evaluate the influence spread of any algorithm using Monte-Carlo simulation, run:

		./IMPC -data NetHEPT.txt -seeds NetHEPT_seeds.txt -k 50
	
	Arguments:
		-data:
			the graph file
		-seeds:
			the seed file contains the ids of seed nodes, separated by space.

- Graph file format

	The first line indicates the number of nodes and edges.
	The following lines includes the source and destination node ids of an edge, together with the propagation probability

	line 1 : #nodes	#edges
	line 2 to 1+#edge : src_id    dest_id    propagation_probability

	All inputs are separated by tab(\t).

	For undirected graph, please keep two links in double direction.

	Example:
	4	4
	0	1	0.02
	1	0	0.05
	2	3	0.1
	3	2	0.2

	This graph contains four nodes and four edges. 

	Sample graph file "NetHEPT.txt" is included.


- Community file format

	The file contains |C| lines, where |C| is the number of communities. 
	Each line represents a community, consisting of ids (separated by \t) of the nodes in that community

	Example:
	0	1
	2	3

	The network consists of two communities {0,1} and {2,3}.

	Sample community file "NetHEPT_com.txt" is included.

To run the program correctly, please make sure that all node ids are ranged from 0 to #nodes-1.
Also make sure that the ids in the community file is consist with the ids in the graph file.
