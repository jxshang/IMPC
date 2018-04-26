#include "Utility.h"
#include "Graph.h"
#include "MemoryUsage.h"

#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<string>
#include<vector>
#include<map>
#include<set>
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
void run_impc(Graph *g, string data, string com, int k, float alpha, float beta);
double *compute_potential(Graph *g);
double get_com_inf(int num_nodes, int com_size, float alpha, float beta);
double get_inf(Graph *g,  int node, double *ptl_arr, int *n2c, int *c2size, float alpha, float beta);
double marginal_gain(Graph *g, set<int> seed_set, int node, double total_score, double *ptl_arr, double *r_arr,
					 set<int> neigh_set, int *n2c, int *c2size, map<int, int> nc_map, float alpha, float beta);
void add_seed(Graph *g, set<int> seed_set, int node, double *r_arr, set<int> neigh_set, int *n2c, map<int, int> nc_map);

void evaluate(Graph *g, string data, string seeds, int k);
void evaluate_total(Graph *g, string data, string seeds, int k, int simus);
float mc_influence(Graph *g, int *seed_arr, int k);
float mc_influence(Graph *g, int *seed_arr, int k, int simus);

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
    string data = "";  // the path of the dataset
	string data_out = ""; //used to convert graph type
    int k=0;  //the # of seeds to be found
    string com = "";
	float alpha = -1;   //the algorithm parameter
	float beta = -1;
    string seeds = "";  // the path of the seed nodes for MC simulation
	bool is_total = false, is_stat = false;
	int simus = 10000;
	int diff_param_type = GIVEN_PARAM;  // the diffusion model: WC, IC, and TRIVALENCY
	vector<double> diff_param;  // the diffusion model parameters

    for(int i=0; i<argn; i++)
    {
        if(argv[i] == string("-data"))
        	data = string(argv[i+1]);
        if(argv[i] == string("-k"))
        	k = atoi(argv[i+1]);
        if(argv[i] == string("-com"))
        	com = argv[i+1];
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
		if(argv[i] == string("-WC"))
			diff_param_type = WC_PARAM;
		if(argv[i] == string("-pp")){ // e.g., -pp 0.1,0.01,0.001
			vector<string> str_vec;
			split(str_vec, string(argv[i+1]), ",");
			for(int i = 0; i < str_vec.size(); i++){
				double pp = atof(str_vec[i].c_str());
				diff_param.push_back(pp);
			}
			if(diff_param.size() == 1)
				diff_param_type = IC_PARAM;
			else
				diff_param_type = TRIVAL_PARAM;
		}
    }
    if (data=="")
        ExitMessage("argument data missing");
    Graph *g = new Graph(data);
	g->genProbT(diff_param_type, diff_param);
	cout << "graph " << data << " was built successfully!" << endl;
	disp_mem_usage("");
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
	if (seeds == ""){
		run_impc(g, data, com, k, alpha, beta);
	}
    else{
		if(is_total)
			evaluate_total(g, data, seeds, k, simus);
		else
			evaluate(g, data, seeds, k);
	}
}

void run_impc(Graph *g, string data, string com, int k, float alpha, float beta){
	cout << "--------------------------------------------------------------------------------" << endl;
    cout << "data:" << data << " com:" << com << " k:" << k;
	cout << " alpha:" << alpha << " beta:" << beta << endl;
	ifstream ifs(com.c_str());
	if (!ifs){
		cout << "community file: " << com << " not openned!" << endl;
		return;
	}	
	else
		cout << "community file: " << com << " successfully opened!" << endl;
	int *n2c = new int[g->n];
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
	for(int i = 0; i < g->n; i++){
		c2size[n2c[i]]++;
	}
	cout << "# of communities: " << com_id << endl;
	clock_t time_start = clock();
	cout << "Finding top " << k << " nodes with IMPC algorithm" << endl;
	cout << "No.\tnode_id\ttime(s)\tScore" << endl;

	//Compute optential for all nodes
	double *ptl_arr = compute_potential(g);

	//Initialize the reverse array
	double *r_arr = new double[g->n];
	for(int i = 0; i < g->n; i++){
		r_arr[i] = 1.0;
	}

	int *seed_arr = new int[k];
	float *score_arr = new float[k];
	double *time_arr = new double[k];

	priority_queue<MinPair, vector<MinPair>, greater<MinPair> > tmp_pqueue;

	set<int> tmp_set;
	for (int i = 0; i < g->n; i++){
		tmp_set.insert(i);
		double score = get_inf(g, i, ptl_arr, n2c, c2size, alpha, beta);
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
	set<int> neigh_set;  // the neighbor set of S, i.e., N(S)
	map<int, int> nc_map; // the neighbor community map
	
	int *updated = new int[g->n];  //the flag array indicates whehter marginal gain of a node is updated
	for (int i = 0; i < g->n; i++)
		updated[i] = 1;
	double total_score = 0;
	for (int i = 0; i < k; i++) {
		Pair best_pair = pqueue.top();
		pqueue.pop();
		//Find the best candidate with CELF strategy
		while (!updated[best_pair.key]) {
			double increase = marginal_gain(g, seed_set, best_pair.key, total_score, ptl_arr, r_arr, neigh_set, n2c, c2size, nc_map, alpha, beta);
			best_pair.value = increase;
			updated[best_pair.key] = 1;
			pqueue.push(best_pair);
			best_pair = pqueue.top();
			pqueue.pop();
		}
		add_seed(g, seed_set, best_pair.key, r_arr, neigh_set, n2c, nc_map);
		total_score += best_pair.value;
		seed_arr[i] = best_pair.key;
		score_arr[i] = total_score;
		time_arr[i] = (double)(clock() - time_start) / CLOCKS_PER_SEC;
		cout << i+1 << "\t" << best_pair.key << "\t";
		printf("%.4f\t%.6f\n", time_arr[i], total_score);
		memset(updated, 0, g->n * sizeof(int)); // reset the flag array
	}

	cout <<"Scores:";
	for(int i = 0; i < k; i++){
		cout << " " << score_arr[i];
	}
	cout << endl;

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

double *compute_potential(Graph *g){
	double *ptl_arr = new double[g->n];
	for(int i = 0; i < g->n; i++){
		ptl_arr[i] = 0;
		int k_out = g->gT[i].size();
		for(int j = 0; j < k_out; j++){
			ptl_arr[i] += g->probT[i][j];
		}
	}
	return ptl_arr;
}

double get_com_inf(int num_nodes, int com_size, float alpha, float beta){
	return alpha * (1-beta*exp(-2.0*com_size*num_nodes / (com_size+num_nodes)));
}

double get_inf(Graph *g, int node, double *ptl_arr, int *n2c, int *c2size, float alpha, float beta){
	double inf = 1; //the influence of node on itself
	inf += ptl_arr[node];  //the influence of node on its direct neighbors
	int k_out = g->gT[node].size();
	map<int, int> nc_map;
	for(int i = 0; i < k_out; i++){
		int neigh = g->gT[node][i];
		double inf_prob = g->probT[node][i];
		inf += inf_prob * ptl_arr[neigh];
		int neigh_com = n2c[neigh];
		if(nc_map.find(neigh_com) == nc_map.end())
			nc_map.insert(pair<int, int>(neigh_com, 1));
		else
			nc_map[neigh_com] = nc_map[neigh_com]+1;
	}
	// compute the influence within communities
	map<int, int>::iterator it;
	for(it = nc_map.begin(); it != nc_map.end(); it++){
		int neigh_com = it->first;
		int num_nodes = it->second;
		int com_size = c2size[neigh_com];
		inf += get_com_inf(num_nodes, com_size, alpha, beta);
	}
	return inf;
}

double marginal_gain(Graph *g, set<int> seed_set, int node, double total_score, double *ptl_arr, double *r_arr, set<int> neigh_set, int *n2c, int *c2size, map<int, int> nc_map, float alpha, float beta){
	double gain = 1.0; //the increase by node itself
	gain += ptl_arr[node]; //the increase by node's potential
	if(neigh_set.find(node) != neigh_set.end()){ // if node is in the neighbor set of S
		gain -= (1-r_arr[node])*ptl_arr[node]; //its previously computed influence should be eliminated
	}
	int k_out = g->gT[node].size();
	map<int, int> tmp_map;
	for(int i = 0; i < k_out; i++){
		int neigh = g->gT[node][i]; // the neigh of node
		gain += r_arr[neigh] * g->probT[node][i] * ptl_arr[neigh];
		if(neigh_set.find(neigh) == neigh_set.end()){  //if the neigh is new
			int neigh_com = n2c[neigh];
			if(tmp_map.find(neigh_com) == tmp_map.end())
				tmp_map.insert(pair<int, int>(neigh_com, 1));
			else
				tmp_map[neigh_com] = tmp_map[neigh_com] + 1;
		}
	}
	// Compute the marginal gain of community influence
	for(map<int, int>::iterator it = tmp_map.begin(); it != tmp_map.end(); it++){
		int neigh_com = it->first;
		int new_nodes = it->second;
		int old_nodes = nc_map[neigh_com];
		int com_size = c2size[neigh_com];
		gain += get_com_inf(old_nodes+new_nodes, com_size, alpha, beta) - get_com_inf(old_nodes, com_size, alpha, beta);
	}
	return gain;
}

void add_seed(Graph *g, set<int> seed_set, int node, double *r_arr, set<int> neigh_set, int *n2c, map<int, int> nc_map){
	seed_set.insert(node); 
	r_arr[node] = 0.0;
	int k_out = g->gT[node].size();
	for(int i = 0; i < k_out; i++){
		int neigh = g->gT[node][i];
		if(r_arr[neigh] == 1.0){  // if the not influence probability of neigh is 1, then neigh must not in N(S)
			neigh_set.insert(neigh);  //add neigh to N(S)
			int neigh_com = n2c[neigh];  //add neigh_com to NC(S)
			if(nc_map.find(neigh_com) == nc_map.end())
				nc_map.insert(pair<int, int>(neigh_com, 1));
			else
				nc_map[neigh_com] = nc_map[neigh_com]+1;
		}
		double inf_prob = g->probT[node][i]; // the influence probability from node to neigh
		r_arr[neigh] *= (1-inf_prob); //update the not influence probability from S to neigh
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

float mc_influence(Graph *g, int *seed_arr, int k){
	srand((unsigned)time(NULL));
	double inf = 0;
	int *i_arr = new int[g->n]; //the array of current active nodes
	int i_size = 0; // the # of newly active nodes 
	int *r_arr = new int[g->n]; // the array of previous active nodes
	int r_size = 0; // the # of previously active nodes
	int *si_arr = new int[g->n];  // the array of nodes to be active in t+1
	int si_size = 0; // the # of nodes to be active in t+1
	int *state_arr = new int[g->n]; // the state of nodes
	memset(state_arr, S_STATE, g->n * sizeof(int)); // initialize the state array	
	int *rand_arr = new int[g->n]; //the 0 ~ n-1 numbers sorted by random order
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
				int k_out = g->gT[i_node].size();
				for(int j = 0; j < k_out; j++){
					int neigh = g->gT[i_node][j];
					if (state_arr[neigh] == S_STATE) {
						double pp = g->probT[i_node][j];
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

float mc_influence(Graph *g, int *seed_arr, int k, int simus){
	srand((unsigned)time(NULL));
	double inf = 0;
	int *i_arr = new int[g->n]; //the array of current active nodes
	int i_size = 0; // the # of newly active nodes 
	int *r_arr = new int[g->n]; // the array of previous active nodes
	int r_size = 0; // the # of previously active nodes
	int *si_arr = new int[g->n];  // the array of nodes to be active in t+1
	int si_size = 0; // the # of nodes to be active in t+1
	int *state_arr = new int[g->n]; // the state of nodes
	memset(state_arr, S_STATE, g->n * sizeof(int)); // initialize the state array	
	int *rand_arr = new int[g->n]; //the 0 ~ n-1 numbers sorted by random order
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
				int k_out = g->gT[i_node].size();
				for(int j = 0; j < k_out; j++){
					int neigh = g->gT[i_node][j];
					if (state_arr[neigh] == S_STATE) {
						double pp = g->probT[i_node][j];
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
