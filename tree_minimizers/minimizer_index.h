#ifndef MST_SET_H
#define MST_SET_H

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <list>
#include "heavy_string.h"


using namespace std;

class MinimizerIndex{
	map<char, int> amap;
	int N;
	string alph;	
	vector<vector<double>> fP;

	string H;
	vector<double> pi_prefix;
	
	friend std::istream & operator >> (std::istream& input, MinimizerIndex &M);
	int pruning(int first_pos,double p1);
	
public:
	void build_index(double z,int ell);
};

#endif //MST_SET_H  
