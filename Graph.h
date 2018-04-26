#pragma once
#include<string>
#include<iostream>
#include<fstream>
#include<vector>

#define GIVEN_PARAM 0
#define WC_PARAM 1
#define IC_PARAM 2
#define TRIVAL_PARAM 3

using namespace std;

class Graph{
public:
	int n;
	int m;
	vector<vector<int> > gT;
	vector<vector<double> > probT;
	
	Graph(string graph_file){
		cout << "Building graph from rom: " << graph_file << " ..." << endl;
		ifstream ifs(graph_file.c_str());

		if(!ifs)
			cout << "graph file: " << graph_file << " not openned!" << endl;
		else
			cout << "graph file: " << graph_file << " successfully opened!" << endl;

		string buffer1, buffer2, buffer3;

		ifs >> buffer1 >> buffer2;
		n = atoi(buffer1.c_str());
		m = atoi(buffer2.c_str());
		cout << "(#nodes, #edges): (" << n << ", " << m << ")" << endl;
		//initialize the adjacency matrix
		for(int i = 0; i < n; i++){
			gT.push_back(vector<int>());
			probT.push_back(vector<double>());
		}

		cout << "reading edges ..." << endl;
		for(int i = 0; i < m; i++){
			ifs >> buffer1 >> buffer2 >> buffer3;
			int src = atoi(buffer1.c_str());
			int dest = atoi(buffer2.c_str());
			double prob = atof(buffer3.c_str());

			gT[src].push_back(dest);
			probT[src].push_back(prob);
		}
		ifs.close();
	}

	void print(){
		cout << "print graph with " << n << " nodes and " << m << " edges" << endl;
		for(int i = 0; i < m; i++){
			for(int j = 0; j < gT[i].size(); j++){	
				cout << i << " "  << gT[i][j] << " " << probT[i][j] << endl;
			}
		}
	}

	void genProbT(int diff_param_type, vector<double> params){
		if(diff_param_type == WC_PARAM){
			vector<int> inDeg;
			for(int i = 0; i < n; i++){
				inDeg.push_back(gT[i].size());
			}
			for(int i = 0; i < n; i++){
				int neighbors = gT[i].size();
				for(int j = 0; j < neighbors; j++){
					int dest = gT[i][j];
					probT[i][j] = 1.0 / inDeg[dest];
				}
			}
			inDeg.clear();
		}
		else if(diff_param_type == IC_PARAM){
			double prob = params[0];
			for(int i = 0; i < n; i++){
				int neighbors = gT[i].size();
				for(int j = 0; j < neighbors; j++){
					probT[i][j] = prob;
				}
			}
		}
		else if(diff_param_type == TRIVAL_PARAM){
			srand(time(NULL));
			int num = params.size();
			for(int i = 0; i < n; i++){
				int neighbors = gT[i].size();
				for(int j = 0; j < neighbors; j++){
					double prob = params[rand()% num];
					probT[i][j] = prob;
				}
			}
		}
	}

	void write_graph(string graph_out){
		ofstream ofs(graph_out.c_str());
		ofs << n << " " << m << endl;
		for(int i = 0; i < n; i++){
			int neighbors = gT[i].size();
			for(int j = 0; j < neighbors; j++){
				ofs << i << "\t" << gT[i][j] << "\t" << probT[i][j] << endl;
			}
		}
		ofs.close();
	}
};
