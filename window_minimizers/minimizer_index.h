#pragma once

#include<string>
#include<vector>
#include<iostream>

class MinimizerIndex {
	int N;
	int Nz;
	std::string alph;
	std::vector<std::vector<double>> fP;

	friend std::istream & operator >> (std::istream& input, MinimizerIndex &M);

public:
	MinimizerIndex(): alph(), fP(){}
	void build_index(double z, int ell);
};
