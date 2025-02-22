#include <cmath>
#include <chrono>
#include <list>
#include <cstdlib>
#include <ctime>
#include <set>
#include <deque>
#include <unordered_map>
#include <sys/resource.h>
#include <algorithm>
#include <unordered_map>

#include "estimation.h"
#include "krfp.h"
#include "minimizers.h"
#include "minimizer_queues.h"
#include "minimizer_index.h"

using std::endl;
using std::cerr;
using get_time = std::chrono::steady_clock;


bool sort_sa(const pair<int,int> &a,const pair<int,int> &b)
{
       return a.first<b.first;
}

std::istream & operator >> (std::istream& input, MinimizerIndex &M) {
    input >> M.N;
    input >> M.alph;
    int A = M.alph.size();
    for (int i = 0; i < M.N; ++i) {
        double sum = 0;
        std::vector<double> symbol(A, 0);
        for (int j = 0; j < A; ++j) {
            input >> symbol[j];
            sum += symbol[j];
        }
        if (std::abs(sum-1) > EPS) {
            std::cerr << "Probabilities at position " << i << " do not sum up to 1" << std::endl;
            throw 1;
        }
        M.fP.emplace_back(symbol);
    }
 	return input;
}

std::chrono::duration<double> test_minimizer_queue(MinimizerQueue &mq, string &text, INT ell){
	INT n= text.length();
	auto begin = get_time::now();
	vector<INT> minimizers={0};
	for(INT j = 0; j < ell - 1; ++j){
		mq.add_right(text[j]);
	}
	for(INT j = ell - 1; j < n; ++j){
		mq.add_right(text[j]);
		if(minimizers.back() != mq.min_pos())
			minimizers.push_back(mq.min_pos());		
		mq.rem_left();
	}
	//cout << "Different minimizers count: " << minimizers.size()<< endl;
	auto end = get_time::now();
	return end - begin;
}


void MinimizerIndex::build_index(double z, int ell){
	vector<vector<double>> rP(fP.rbegin(), fP.rend());
	Estimation fS(fP,alph,z);
	PropertyString fT;
	std::vector<int> f_mini_pos;
	std::vector<int> r_mini_pos;
	int i = 0;
	int k = ceil(4 * log2(ell) / log2(alph.size()));
	karp_rabin_hashing::init(ceil(4 * log2(ell) / log2(alph.size())));
	int n = fS.strings()[0].length();
	cout << "k " << k << endl;
	int w = ell - k + 1;
	
	cout << "Operations count (nz): " << n * ((INT) z)<<endl;	
	
	auto begin = get_time::now();
	for(PropertyString const & s : fS.strings()){
		std::unordered_set<uint64_t> M;
		string temp_s = s.string();
		compute_minimizers(temp_s, w, k, M);
		i++;
	}
	auto end = get_time::now();
	auto diff2 = end - begin;
	cout << "Sliding window (original) time: "<< chrono::duration_cast<chrono::milliseconds>(diff2).count()<<endl;



	begin = get_time::now();
	for(PropertyString const & s : fS.strings()){
		string temp_s = s.string();
		compute_minimizers2(temp_s, w, k);
		i++;
	}

	end = get_time::now();
	diff2 = end - begin;
	cout << "Sliding window (vector) time: "<< chrono::duration_cast<chrono::milliseconds>(diff2).count()<<endl;



	HeapMinimizer hm(k);
	TwoStackMinimizer tsm(k);	
	std::chrono::duration<double> time1;
	for(PropertyString const & s : fS.strings()){
		string temp_s = s.string();
		time1 += test_minimizer_queue(hm, temp_s, ell);
	}
	cout << "Heap minimizer time: " << chrono::duration_cast<chrono::milliseconds>(time1).count() << endl;

	std::chrono::duration<double> time2;
	for(PropertyString const & s : fS.strings()){
		string temp_s = s.string();
		time2 += test_minimizer_queue(tsm, temp_s, ell);
	}
	cout << "Two stack minimizer time: " << chrono::duration_cast<chrono::milliseconds>(time2).count() << endl;	
	
}
