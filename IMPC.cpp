#include "Utility.h"
#include "Graph.h"
#include "MemoryUsage.h"

#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<string>
#include<vector>
#include<map>
#include<fstream>
#include<sstream>
#include<queue>
#include<functional>

#define S_STATE 0
#define I_STATE 1
#define SI_STATE 2
#define R_STATE 3
#define REDUND 10

struct Pair {
	int key;
	double value;
	Pair(int key, double value) :key(key), value(value) {};
};
typedef struct Pair Pair;

bool operator > (Pair a, Pair b) {
	return a.value < b.value;
}

struct MinPair {
	int key;
	double value;
	MinPair(int key, double value) :key(key), value(value) {};
};
typedef struct MinPair MinPair;

bool operator > (MinPair a, MinPair b) {
	return a.value > b.value;
}

using namespace std;

void parseArg(int argn, char ** argv);
void run_impc(Graph *g, string data, string com, int k, float alpha, float beta, int level);
void evaluate(Graph *g, string data, string seeds, int k);
void evaluate_total(Graph *g, string data, string seeds, int k, int simus);

float mc_influence(Graph *g, int *seed_arr, int k);
float mc_influence(Graph *g, int *seed_arr, int k, int simus);
double get_score(Graph *g, set<int> seed_set, int *n2c, int *c2size, float alpha, float beta);
double get_score_1(Graph *g, set<int> seed_set, int *n2c, int *c2size, float alpha, float beta);

int main(int argn, char ** argv)
{
	cout << "Program Start at: " << currentTimestampStr() << endl;
	cout << "Arguments: ";
	for(int i = 0; i < argn; i++){
		cout << argv[i]<<" ";
	}
	cout << endl;
	cout << "--------------------------------------------------------------------------------" << endl;
    parseArg( argn, argv );
    cout << "--------------------------------------------------------------------------------" << endl;
    cout<<"Program Ends Successfully at: " << currentTimestampStr() << endl;
    return 0;
}

void parseArg(int argn, char ** argv)
{
	// the parameters
    string data="";  // the path of the dataset
    int k=0;  //the # of seeds to be found
    string com = "";
	int level = 2;
	float alpha = -1;   //the algorithm parameter
	float beta = -1;
    string seeds = "";  // the path of the seed nodes for MC simulation
	bool is_total = false, is_stat = false;
	int simus = 10000;

    for(int i=0; i<argn; i++)
    {
        if(argv[i] == string("-data"))
        	data = string(argv[i+1]);
        if(argv[i] == string("-k"))
        	k = atoi(argv[i+1]);
        if(argv[i] == string("-com"))
        	com = argv[i+1];
		if(argv[i] == string("-level"))
			level = atoi(argv[i+1]);
		if(argv[i] == string("-alpha"))
			alpha = atof(argv[i+1]);
		if(argv[i] == string("-beta"))
			beta = atof(argv[i+1]);
        if(argv[i] == string("-seeds"))
        	seeds = argv[i+1];
		if(argv[i] == string("-total"))
			is_total = true;
		if(argv[i] == string("-simus"))
			simus = atoi(argv[i+1]);
    }
    if (data=="")
        ExitMessage("argument data missing");
	if(k == 0 && !is_stat)
		ExitMessage("argument k missing");
    if(seeds == "" && !is_stat){
    	if(com == "")
    		ExitMessage("argument com is missing");
		if(alpha == -1)
			ExitMessage("argument alpha is missing");
		if(beta == -1)
			ExitMessage("argument beta is missing");
    }	
    Graph *g = new Graph(data);
	cout << "graph " << data << " was built successfully!" << endl;
	if (seeds == ""){
		run_impc(g, data, com, k, alpha, beta, level);
	}
    else{
		if(is_total)
			evaluate_total(g, data, seeds, k, simus);
		else
			evaluate(g, data, seeds, k);
	}
}

void evaluate(Graph *g, string data, string seeds, int k){
	cout << "evaluating influence... data:" << data << " seeds:" << seeds << endl;
	int *seed_arr = new int[k];
	ifstream ifs(seeds.c_str());
	if(!ifs)
		cout << "seeds file: " << seeds << " cannot be openned!" << endl;
	else
		cout << "seeds file: " << seeds << " successfully openned!" << endl;
	cout << "id\tseed\tinfluence\ttimestamp" << endl;
	string buffer;
	int point_arr[11] = {1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
	float *inf_arr = new float[11];
	int id = 0;
	for(int i = 0; i < k; i++){
		ifs >> buffer;
		seed_arr[i] = atoi(buffer.c_str());
		int match = 0;
		for(int j = 0; j < 11; j++){
			if(point_arr[j] == i+1){
				match = 1;
				break;
			}
		}
		if(match){
			float inf = mc_influence(g, seed_arr, i + 1);
			inf_arr[id++] = inf;
			cout << i + 1 << "\t" << seed_arr[i] << "\t";
			cout << inf << '\t' << currentTimestampStr() << endl;
		}
	}
	cout << "inf=[" << inf_arr[0];
	for(int i = 1; i < id; i++)
		cout << ", " << inf_arr[i];
	cout << "];" << endl;
}

float mc_influence(Graph *g, int *seed_arr, int k){
	srand((unsigned)time(NULL));
	double inf = 0;
	int *i_arr = new int[g->num_nodes]; //the array of current active nodes
	int i_size = 0; // the # of newly active nodes 
	int *r_arr = new int[g->num_nodes]; // the array of previous active nodes
	int r_size = 0; // the # of previously active nodes
	int *si_arr = new int[g->num_nodes];  // the array of nodes to be active in t+1
	int si_size = 0; // the # of nodes to be active in t+1
	int *state_arr = new int[g->num_nodes]; // the state of nodes
	memset(state_arr, S_STATE, g->num_nodes * sizeof(int)); // initialize the state array	
	int *rand_arr = new int[g->num_nodes]; //the 0 ~ n-1 numbers sorted by random order
	for(int r = 0; r < NUM_SIMUS; r++){
		double active_size = 0;
		//reset the state of all nodes		
		for(int i = 0; i < r_size; i++){
			state_arr[r_arr[i]] = S_STATE;
		}		
		r_size = 0;		
		// initialize the seed nodes
		for(int i = 0; i < k; i++){
			i_arr[i_size++] = seed_arr[i];
			state_arr[i_arr[i]] = I_STATE;
		}
		while(i_size > 0){
			active_size += i_size;
			si_size = 0;
			randomOrder(rand_arr, i_size);
			for(int i = 0; i < i_size; i++){
				int i_node = i_arr[i];
				int k_out = g->node_array[i_node].k_out;
				for(int j = 0; j < k_out; j++){
					int neigh = g->node_array[i_node].id_array[j];
					if (state_arr[neigh] == S_STATE) {
						int k_in = g->node_array[neigh].k_in;
						double pp = 1.0 / k_in;
						double rand_float = ((double)rand()) / RAND_MAX;
						if(rand_float < pp) {
							state_arr[neigh] = SI_STATE;
							si_arr[si_size++] = neigh;
						}
					}					
				}
			}
			for(int i = 0; i < i_size; i++){
				state_arr[i_arr[i]] = R_STATE;
				r_arr[r_size++] = i_arr[i];
			}
			i_size = 0;
			for(int i = 0; i < si_size; i++){
				state_arr[si_arr[i]] = I_STATE;
				i_arr[i_size++] = si_arr[i];
			}
		}
		inf += active_size;
	}
	delete[] i_arr;
	delete[] r_arr;
	delete[] si_arr;
	delete[] state_arr;
	delete[] rand_arr;
	return inf / NUM_SIMUS;
}

void run_impc(Graph *g, string data, string com, int k, float alpha, float beta, int level){
	cout << "--------------------------------------------------------------------------------" << endl;
	cout << "Start IMPC algorithm" << endl;
    cout << "data:" << data << " com:" << com << " k:" << k;
	cout << " alpha:" << alpha << " beta:" << beta << endl;
	ifstream ifs(com.c_str());
	if (!ifs){
		cout << "community file: " << com << " not openned!" << endl;
		return;
	}	
	else
		cout << "community file: " << com << " successfully opened!" << endl;
	int *n2c = new int[g->num_nodes];
	string str;
	int com_id = 0;
	while (getline(ifs, str)){
		istringstream iss(str);
		string buffer;
		while (iss >> buffer) {
			int node = atoi(buffer.c_str());
			n2c[node] = com_id;
		}
		com_id++;
	}
	int* c2size = new int[com_id];
	memset(c2size, 0, com_id*sizeof(int));
	for(int i = 0; i < g->num_nodes; i++){
		c2size[n2c[i]]++;
	}
	cout << "# of communities: " << com_id << endl;
	clock_t time_start = clock_t();
	cout << "Finding top " << k << " nodes with IMPC algorithm" << endl;
	cout << "No.\tnode_id\ttime(s)\tScore" << endl;

	int *seed_arr = new int[k];
	float *score_arr = new float[k];
	double *time_arr = new double[k];

	priority_queue<MinPair, vector<MinPair>, greater<MinPair> > tmp_pqueue;

	set<int> tmp_set;
	for (int i = 0; i < g->num_nodes; i++){
		tmp_set.insert(i);
		double score = 0;
		if(level == 1)
			score = get_score_1(g, tmp_set, n2c, c2size, alpha, beta);
		else
			score = get_score(g, tmp_set, n2c, c2size, alpha, beta);
		tmp_set.erase(tmp_set.begin());
		MinPair m_pair(i, score);
		if ((int)tmp_pqueue.size() >= REDUND * k && score <= tmp_pqueue.top().value)
			continue;
		tmp_pqueue.push(m_pair);
		if ((int)tmp_pqueue.size() > REDUND * k) {  // Kepp the top 10*k nodes with the maximum score
			tmp_pqueue.pop();
		}
	}

	priority_queue<Pair, vector<Pair>, greater<Pair> > pqueue;
	while (!tmp_pqueue.empty()) {
		MinPair min_pair = tmp_pqueue.top();
		tmp_pqueue.pop();
		Pair pair(min_pair.key, min_pair.value);
		pqueue.push(pair);
	}

	// find the top k nodes
	set<int> seed_set;  //seed set to be found
	
	int *updated = new int[g->num_nodes];  //the flag array indicates whehter marginal gain of a node is updated
	for (int i = 0; i < g->num_nodes; i++)
		updated[i] = 1;
	double total_score = 0;
	for (int i = 0; i < k; i++) {
		Pair best_pair = pqueue.top();
		pqueue.pop();
		//Find the best candidate with CELF strategy
		while (!updated[best_pair.key]) {
			seed_set.insert(best_pair.key);
			double new_score = 0;
			if(level == 1)
				new_score = get_score_1(g, seed_set, n2c, c2size, alpha, beta);
			else
				new_score = get_score(g, seed_set, n2c, c2size, alpha, beta);
			seed_set.erase(best_pair.key);
			double increase = new_score - total_score;
			best_pair.value = increase;
			updated[best_pair.key] = 1;
			pqueue.push(best_pair);
			best_pair = pqueue.top();
			pqueue.pop();
		}
		seed_set.insert(best_pair.key);
		total_score += best_pair.value;
		seed_arr[i] = best_pair.key;
		score_arr[i] = total_score;
		time_arr[i] = (double)(clock() - time_start) / CLOCKS_PER_SEC;
		cout << i+1 << "\t" << best_pair.key << "\t";
		printf("%.4f\t%.6f\n", time_arr[i], total_score);
		memset(updated, 0, g->num_nodes * sizeof(int)); // reset the flag array
	}

	cout <<"Seeds:";
	for(int i = 0; i < k; i++){
		cout << " " << seed_arr[i];
	}
	cout << endl;
	
	delete[] n2c;
	delete[] updated;
	delete[] c2size;
	delete[] seed_arr;
	delete[] score_arr;
	delete[] time_arr;
	disp_mem_usage("");
	cout << "Time used: " << (double)(clock() - time_start) / CLOCKS_PER_SEC << " s" << endl;
}

double get_score_1(Graph *g, set<int> seed_set, int *n2c, int *c2size, float alpha, float beta){
	double score = 0;
	set<int>::iterator it;
	set<int> com_set;
	map<int, set<int> > s1addcom;
	set<int> set1;
	map<int, double> mp1;

	double potential = 0;

	// the first layer
	for(it = seed_set.begin(); it != seed_set.end(); it++){
		int seed = *it;
		int *p_neigh = g->node_array[seed].id_array;
		int k_out = g->node_array[seed].k_out;
		for (int i = 0; i < k_out; i++){
			int neigh = p_neigh[i];
			if(set1.find(neigh) != set1.end()){
				mp1[neigh] = mp1[neigh] + 1.0/g->node_array[neigh].k_out;
			}
			else{
				if(seed_set.find(neigh) == seed_set.end()){
					set1.insert(neigh);
					mp1.insert(pair<int, double>(neigh, 1.0/g->node_array[neigh].k_out));
				}
			}
		}
	}

	int temp = 0;
	int temp_com = 0;
	for(it = set1.begin(); it != set1.end(); it++){
		temp = *it;
		temp_com = n2c[temp];
		if(com_set.find(temp_com) == com_set.end()){
			com_set.insert(temp_com);
			set<int> temp_set;
			temp_set.insert(temp);
			s1addcom.insert(pair<int, set<int> >(temp_com, temp_set));
		}
		else{
			s1addcom[temp_com].insert(temp);
		}
		potential += mp1[temp];
	}

	// compute the influence of S1 within communities
	for(it = com_set.begin(); it != com_set.end(); it++){
		temp_com = *it;
		score += alpha * (1-beta*exp(-2.0*c2size[temp_com]*s1addcom[temp_com].size() / (c2size[temp_com]+s1addcom[temp_com].size())));
	}
	score += potential + seed_set.size();
	return score;
}

double get_score(Graph *g, set<int> seed_set, int *n2c, int *c2size, float alpha, float beta){
	double score = 0;
	set<int>::iterator it;
	set<int> com_set;
	map<int, set<int> > s2addcom;
	set<int> set1;
	map<int, double> mp1;
	set<int> set2;
	map<int, double> mp2;

	double potential = 0;

	// the first layer
	for(it = seed_set.begin(); it != seed_set.end(); it++){
		int seed = *it;
		int *p_neigh = g->node_array[seed].id_array;
		int k_out = g->node_array[seed].k_out;
		for (int i = 0; i < k_out; i++){
			int neigh = p_neigh[i];
			if(set1.find(neigh) != set1.end()){
				mp1[neigh] = mp1[neigh] + 1.0/g->node_array[neigh].k_out;
			}
			else{
				if(seed_set.find(neigh) == seed_set.end()){
					set1.insert(neigh);
					mp1.insert(pair<int, double>(neigh, 1.0/g->node_array[neigh].k_out));
				}
			}
		}
	}

	// the second layer
	int temp = 0;
	for(it = set1.begin(); it != set1.end(); it++){
		temp = *it;
		int *p_neigh = g->node_array[temp].id_array;
		int k_out = g->node_array[temp].k_out;
		for(int i = 0; i < k_out; i++){
			int neigh = p_neigh[i];
			if(set2.find(neigh) != set2.end()){
				mp2[neigh] = mp2[neigh] + mp1[temp]*1.0/g->node_array[neigh].k_out;
			}
			else{
				if(seed_set.find(neigh) == seed_set.end()){
					set2.insert(neigh);
					mp2.insert(pair<int, double>(neigh, mp1[temp]*1.0/g->node_array[neigh].k_out));
				}
			}
		}
		potential += mp1[temp];
	}

	int temp_com = 0;
	for(it = set2.begin(); it != set2.end(); it++){
		temp = *it;
		temp_com = n2c[temp];
		if(com_set.find(temp_com) == com_set.end()){
			com_set.insert(temp_com);
			set<int> temp_set;
			temp_set.insert(temp);
			s2addcom.insert(pair<int, set<int> >(temp_com, temp_set));
		}
		else{
			s2addcom[temp_com].insert(temp);
		}
		potential += mp2[temp];
	}

	// compute the influence of S2 within communities
	for(it = com_set.begin(); it != com_set.end(); it++){
		temp_com = *it;
		score += alpha * (1-beta*exp(-2.0*c2size[temp_com]*s2addcom[temp_com].size() / (c2size[temp_com]+s2addcom[temp_com].size())));
	}
	score += potential + seed_set.size();
	return score;
}

void evaluate_total(Graph *g, string data, string seeds, int k, int R){
	cout << "evaluating overall influence... data:" << data << " seeds:" << seeds << endl;
	int *seed_arr = new int[k];
	ifstream ifs(seeds.c_str());
	if(!ifs)
		cout << "seeds file: " << seeds << " cannot be openned!" << endl;
	else
		cout << "seeds file: " << seeds << " successfully openned!" << endl;
	string buffer;
	cout << "Seeds:";
	for(int i = 0; i < k; i++){
		ifs >> buffer;
		seed_arr[i] = atoi(buffer.c_str());
		cout << " " << seed_arr[i];
	}
	cout << endl;
	float total_inf = mc_influence(g, seed_arr, k, R);
	cout << "Total influence: " << total_inf << endl;
}

float mc_influence(Graph *g, int *seed_arr, int k, int simus){
	srand((unsigned)time(NULL));
	double inf = 0;
	int *i_arr = new int[g->num_nodes]; //the array of current active nodes
	int i_size = 0; // the # of newly active nodes 
	int *r_arr = new int[g->num_nodes]; // the array of previous active nodes
	int r_size = 0; // the # of previously active nodes
	int *si_arr = new int[g->num_nodes];  // the array of nodes to be active in t+1
	int si_size = 0; // the # of nodes to be active in t+1
	int *state_arr = new int[g->num_nodes]; // the state of nodes
	memset(state_arr, S_STATE, g->num_nodes * sizeof(int)); // initialize the state array	
	int *rand_arr = new int[g->num_nodes]; //the 0 ~ n-1 numbers sorted by random order
	for(int r = 0; r < simus; r++){
		double active_size = 0;
		//reset the state of all nodes		
		for(int i = 0; i < r_size; i++){
			state_arr[r_arr[i]] = S_STATE;
		}		
		r_size = 0;		
		// initialize the seed nodes
		for(int i = 0; i < k; i++){
			i_arr[i_size++] = seed_arr[i];
			state_arr[i_arr[i]] = I_STATE;
		}
		while(i_size > 0){
			active_size += i_size;
			si_size = 0;
			randomOrder(rand_arr, i_size);
			for(int i = 0; i < i_size; i++){
				int i_node = i_arr[i];
				int k_out = g->node_array[i_node].k_out;
				for(int j = 0; j < k_out; j++){
					int neigh = g->node_array[i_node].id_array[j];
					if (state_arr[neigh] == S_STATE) {
						int k_in = g->node_array[neigh].k_in;
						double pp = 1.0 / k_in;
						double rand_float = ((double)rand()) / RAND_MAX;
						if(rand_float < pp) {
							state_arr[neigh] = SI_STATE;
							si_arr[si_size++] = neigh;
						}
					}					
				}
			}
			for(int i = 0; i < i_size; i++){
				state_arr[i_arr[i]] = R_STATE;
				r_arr[r_size++] = i_arr[i];
			}
			i_size = 0;
			for(int i = 0; i < si_size; i++){
				state_arr[si_arr[i]] = I_STATE;
				i_arr[i_size++] = si_arr[i];
			}
		}
		inf += active_size;
	}
	delete[] i_arr;
	delete[] r_arr;
	delete[] si_arr;
	delete[] state_arr;
	delete[] rand_arr;
	return inf / simus;
}
