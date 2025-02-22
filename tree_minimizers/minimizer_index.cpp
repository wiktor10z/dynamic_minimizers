#include <cmath>
#include <algorithm>
#include <queue>
#include <set>
#include <map>
#include <list>
#include <chrono>
#include <ctime>
#include <functional>
#include <cassert>
#include <tuple>
#include "minimizer_index.h"
#include "minimizers.h"
#include "minimizer_queues.h"
#include "krfp.h"
#include "heavy_string.h"

using namespace std;
using get_time = chrono::steady_clock;

bool isEqual(double a, double b) {
	double epsilon = 1e-14;
    return std::abs(a - b) <= epsilon * std::max(std::abs(a), std::abs(b));
}


std::istream & operator >> (std::istream& input, MinimizerIndex &M) {
    input >> M.N;
    input >> M.alph;
    int A = M.alph.size();
	M.pi_prefix = vector<double> (M.N, 1);
    for (int i = 0; i < M.N; ++i) {
        double sum = 0;
        vector<double> symbol(A, 0);
        for (int j = 0; j < A; ++j) {
            input >> symbol[j];
            sum += symbol[j];
        }
 		int which_max = max_element(symbol.begin(), symbol.end()) - symbol.begin();
		M.H+=(M.alph[which_max]);
		double pi = symbol[which_max];
		if(i == 0){
			M.pi_prefix[i] = log2(pi);
		}else{
			M.pi_prefix[i] = M.pi_prefix[i-1] + log2(pi);
		}       
        
        if (!isEqual(sum,1)) {
            cerr << "Probabilities at position " << i << " do not sum up to 1" << std::endl;
            throw 1;
        }
        M.fP.emplace_back(symbol);
    }
 	return input;
}

int MinimizerIndex::pruning(int first_pos, double p1){
	int n = H.size();
	double base;
	if(first_pos == 0){
		base = p1;
	}else{
		base = p1 + pi_prefix[first_pos - 1];
	}
	int beg_pos = first_pos - 1,window = 1, end_pos; //beg_pos - the furthest position about which we know that it belongs to the solid factor
	while((beg_pos + window < n)&&(pi_prefix[beg_pos + window] >= base)){
		beg_pos += window;
		window *= 2;
	}
	if(beg_pos+window >= n){
		end_pos = n;
	}else{
		end_pos = beg_pos + window;
	}//now we know that beg_pos belongs to the factor and that end_pos does not
	while(beg_pos < end_pos - 1){
		if(pi_prefix[(beg_pos + end_pos)/2] - pi_prefix[first_pos] >= p1){
			beg_pos=(beg_pos + end_pos) /2;
		}else{
			end_pos=(beg_pos + end_pos) /2;
		}
	}
	return beg_pos;
	
}


std::chrono::duration<double> test_minimizer_queue(MinimizerQueue &mq, vector<tuple<int, bool, char>> operations){
	auto begin = get_time::now();
	vector<INT> minimizers = {0};
	for(auto iter = operations.begin(); iter != operations.end(); ++iter){
		switch(get<0>(*iter)){
			case 0:
				if(minimizers.back() != mq.min_pos())
					minimizers.push_back(mq.min_pos());
				break;
			case 1:
				if(get<1>(*iter)) mq.add_right(get<2>(*iter));
				else mq.add_left(get<2>(*iter));
				break;
			case 2:
				if(get<1>(*iter)) mq.rem_right();
				else mq.rem_left();
				break;		
			case 3:
				if(get<1>(*iter)) mq.move_right(get<2>(*iter));
				else mq.move_left(get<2>(*iter));
				break;
		}
	}
	cout << "Different minimizers count: " << minimizers.size()<< endl;
	auto end = get_time::now();
	return end - begin;
}


string get_text(vector<char> S){
	reverse(S.begin(), S.end());
	string text(S.begin(),S.end());
	return text;
}


//Compute minimizers for a leaf to root path using sliding window technique
std::chrono::duration<double> compute_minimizers_on_path(vector<char> S, INT w, INT k){
	string text = get_text(S);
	auto begin = get_time::now();
	compute_minimizers2(text, w, k);
	auto end = get_time::now();
	return end - begin;
}


void MinimizerIndex::build_index(double z, int l){
	
	vector<tuple<int, bool, char>> minimizer_operations; // 0-min, 1-add, 2-rem, 3-move, true - right, false - left
	
	int k = ceil(4*log2(l) / log2(alph.size()));
	int n = fP.size();
	int w = l - k + 1;
	karp_rabin_hashing::init(k);//initialization with fixed k
	
	list<pair<size_t,size_t>> minimizer_substrings;
	size_t minimizer_count=0;
	
	MinimizerHeap heap(n,l,k);
	cout << "k " << k << endl;
	//HeapMinimizer dm(k);
	
	for(int i = 0; i < alph.size(); i++){
		amap[alph[i]] = i;
	}
	
	double p = 1.0;
	int a = n-1;
	unordered_set<int> minimizers;
	list<pair<int, char>> diff;
	list<list<pair<int, char>>> global_diff;
	int sig1 = -1;
	int pos1 = n;
	std::chrono::duration<double> time1; //Only for SW naive check
	bool leaf; //Only for SW naive check - marks whether the last move was to the left, and hence move to the right goes out of a leaf
	while( a != n ){
		int sig = sig1 + 1;
		if( a >= 0 && sig != alph.size()){
			if( p != 1 || alph[sig] != H[a]){
				if( p * fP[a][sig] * z < 1 ){
					sig1 = sig;
					continue;
				}else{
					p *= fP[a][sig];
				}
			}else{
				pos1 = a;
			}
			heap.left(alph[sig]);
			leaf = true; //Only for SW naive check
			if(heap.S.size()<=l){ // heap is already modified here
				//dm.add_left(alph[sig]);
				minimizer_operations.push_back(make_tuple(1, false, alph[sig]));
			}else{
				//dm.move_left(alph[sig]);
				minimizer_operations.push_back(make_tuple(3, false, alph[sig]));
			}
			//assert(heap.lefthash == dm.lefthash);
			//assert(heap.righthash == dm.righthash);
			if(H[a] != alph[sig]){
				diff.push_front(make_pair(a, alph[sig]));
			}
			if(heap.S.size() >= l){
				double pi_cum = 1;
				if( pos1 <= 0 ){	
					pi_cum = pow(2,pi_prefix[l]);	//if full length S is the heavy string and reach the beginning, directly use pi prefix
				}else{
					pi_cum = p * pow(2, pi_prefix[a+l-1] - pi_prefix[pos1-1]);
				}
				if(pi_cum * z >= 1){
					//assert(heap.top() == dm.min_pos() + n);
					minimizer_operations.push_back(make_tuple(0, false, 0));
					minimizers.insert(heap.top());
				}
			}
			a = a - 1;
			sig1 = -1;
		}else{
			if(leaf && (heap.S.size() >= l)){ //Only for SW naive check
				//time1 += compute_minimizers_on_path(heap.S, w, k); // sliding window computation of minimizers starting in the node
			}
			a = a + 1;//a is set to the position of the letter that is being removed
			if(minimizers.find(a) != minimizers.end()){
				minimizers.erase(a);
				//finding the longest string starting at this minimizer with weight >=1/z
				int pos3=pruning(pos1,-log2(p)-log2(z))+1;
				minimizer_substrings.push_back(make_pair(minimizer_count+(size_t)a,minimizer_count+(size_t)pos3));
				global_diff.push_back(diff);
				minimizer_count+=(size_t)n;
			}
			// removing the first letter and restoring variables
			if ( !diff.empty()){
				if(diff.front().first == a) {
					diff.pop_front();
				}
			}
			if(!isEqual(p,1)){
				p /= fP[a][amap[heap.S[heap.S.size()-1]]];
				if(p>0.7) p=1.0; //fixing p in case of precision errors (p cannot be between 0.5 and 1, hence p>0.7 means p=1)
			}else{
				pos1=a+1;
			}
			if(heap.S.size() > 0){
				sig1 = amap[heap.S[heap.S.size()-1]];
				heap.right();
				leaf = false; //Only for SW naive check
				if(heap.S.size() >= l){ // heap is already modified here
					//dm.move_right(heap.S[heap.S.size() - l]);
					minimizer_operations.push_back(make_tuple(3, true, heap.S[heap.S.size() - l]));
				}else{
					//dm.rem_left();
					minimizer_operations.push_back(make_tuple(2, false, 0));
				}
			}
		}
	}
	
	// Checking the running time of the two semi-dynamic minimizers with just the minimizer operations.
	cout << "Operations count: " << minimizer_operations.size()<<endl;
	
	//cout << "Naive sliding window time: " << chrono::duration_cast<chrono::milliseconds>(time1).count() << endl;
	HeapMinimizer hm(k);
	TwoStackMinimizer tsm(k);	

	std::chrono::duration<double> time2;
	//for(int i = 0; i < 10; ++i){
		time2 += test_minimizer_queue(hm, minimizer_operations);
	//}
	cout << "Heap minimizer time: " << chrono::duration_cast<chrono::milliseconds>(time2).count() << endl;

	std::chrono::duration<double> time3;
	//for(int i = 0; i < 10;  ++i){
		time3 += test_minimizer_queue(tsm, minimizer_operations);
	//}
	cout << "Two stack minimizer time: " << chrono::duration_cast<chrono::milliseconds>(time3).count() << endl;
}
