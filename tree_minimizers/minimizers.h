#include <string>
#include <vector>
#include <utility>
#include <unordered_set>

using namespace std;

class MinimizerHeap{
	private:
		uint64_t n,l,k,len,lefthash,righthash, letter_k;	
		set<pair<uint64_t,uint64_t>> heap; // (KR hash, position)
	public:
	    vector<char> S;
	    MinimizerHeap(uint64_t n1, uint64_t l1, uint64_t k1);
	   
	    uint64_t top();
	    void left(char a);
	    void right();
};

uint64_t linear_minimizer(vector<char>& S, uint64_t l, uint64_t k);

uint64_t compute_minimizers2(std::string& text, uint64_t w, uint64_t k);
