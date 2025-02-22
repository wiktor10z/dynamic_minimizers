#include <iostream>
#include <deque>
#include <set>
#include <stack>
#include "krfp.h"
#include "minimizer_queues.h"

using namespace std;

typedef uint64_t UINT;
typedef int64_t INT;


// Changing the hash of the string lS to the the one of the string Sr
// Letter l = 0 means no letter on the left (extension by one letter) 
// Letter r = 0 means no letter (removal of a letter) on the left
// We assume len <= k

void hashright(UINT& hash, char l, char r, int len){
	if (r != 0){
		hash = karp_rabin_hashing::concat(hash, r, 1);
	}
	if (l != 0){
		hash = karp_rabin_hashing::fast_subtract(hash, l, len);
	}
}

// Changing the hash of the string Sr to the the one of the string lS
// Letter r = 0 means no letter on the right (extension by one letter) 
// Letter l = 0 means no letter (removal of a letter) on the right
// We assume len <= k
void hashleft(UINT& hash, char l, char r, int len){
	if (l != 0){
		hash = karp_rabin_hashing::fast_concat(l, hash, len);
	}
	if (r != 0){
		hash = karp_rabin_hashing::subtract(hash, r, 0);
		hash = karp_rabin_hashing::leftshift(hash);		 
	}
}

void MinimizerQueue::add_right_update_S(char c){
	S.push_back(c);
	if(S.size() > k){
		hashright(righthash, S[S.size() - k - 1], c, k);				
	}else{
		hashright(righthash, 0, c, S.size() - 1);
		lefthash = righthash; // The string is too short to have multiple factors of length k - lefthash = righthash
	}
}

void MinimizerQueue::rem_right_update_S(){
	if(S.size() > k){
		hashleft(righthash, S[S.size() - k - 1], S[S.size()-1], k);
	}else{
		hashleft(righthash, 0, S[S.size()-1], S.size() - 1);
		lefthash = righthash; // The string is too short to have multiple factors of length k - lefthash = righthash
	}
	S.pop_back();
}


void MinimizerQueue::add_left_update_S(char c){
	S.push_front(c);
	start_pos --;
	if(S.size() > k){
		hashleft(lefthash, c, S[k], k);					
	}else{
		hashleft(lefthash, c, 0, S.size() - 1);
		righthash = lefthash; // The string is too short to have a factor of length k - lefthash = righthash
	}
}

void MinimizerQueue::rem_left_update_S(){
	if(S.size() > k){
		hashright(lefthash, S[0], S[k], k);
	}else{
		hashright(lefthash, S[0], 0, S.size() - 1);
		righthash = lefthash; // The string is too short to have multiple factors of length k - lefthash = righthash
	}
	S.pop_front();
	start_pos ++;
}

	
void HeapMinimizer::add_right(char c){
	add_right_update_S(c);
	if(S.size() >= k){
		heap.insert(make_pair(righthash, start_pos + S.size() - k));	
	}
}	
	
void HeapMinimizer::rem_right() {
	if(heap.size() > 0){
		heap.erase(heap.find(make_pair(righthash, start_pos + S.size() - k)));
	}
	rem_right_update_S();
}		
	
void HeapMinimizer::add_left(char c){
	add_left_update_S(c);
	if(S.size() >= k){
		heap.insert(make_pair(lefthash, start_pos));
	}
}

void HeapMinimizer::rem_left(){
	if(heap.size() > 0){
		heap.erase(heap.find(make_pair(lefthash, start_pos)));
	}
	rem_left_update_S();
}
	
INT HeapMinimizer::min_pos(){
	return heap.begin()->second;
}

void TwoStackMinimizer::clear_stacks(){
	stacks[0] = stack<pair<UINT,INT>>(); // stack does not have function clear()
	stacks[1] = stack<pair<UINT,INT>>();	
}

void TwoStackMinimizer::rebuild(){
	clear_stacks();
	UINT hash = lefthash;
	stack<pair<UINT,INT>> temp_left;
	// Reconstruction with separate stack using Theta(n) space
	// We can reduce the space to the actual size of the structure,
	// but at a cost of few extra operations
	for(INT i = 0; i < (S.size() - k + 1) / 2; ++i){
		temp_left.push(make_pair(hash, start_pos + i));
		hashright(hash, S[i], S[i + k], k);
	}
	stacks[0].push(temp_left.top());
	temp_left.pop();
	while(temp_left.size() > 0){
		if(temp_left.top().first <= stacks[0].top().first){
			stacks[0].push(temp_left.top());
		}
		temp_left.pop();
	}
	// Left stack is rebuilt, now we rebuild the right stack
	stacks[1].push(make_pair(hash, start_pos + (S.size() - k + 1) / 2));
	for(INT i = (S.size() - k + 1) / 2; i < S.size() - k ; ++i){
		hashright(hash, S[i], S[i + k], k);	
		if(hash < stacks[1].top().first){
			stacks[1].push(make_pair(hash, start_pos + i + 1));
		}	
	}
}
 
void TwoStackMinimizer::add_right(char c){
	add_right_update_S(c);
	if(S.size() >= k){
		if((stacks[1].size() == 0) || (stacks[1].top().first > righthash)){
			stacks[1].push(make_pair(righthash,start_pos + S.size() - k));
		}
	}
}	


void TwoStackMinimizer::rem_right() {
	if(S.size() > k){
		if(stacks[1].size() == 0){ // the right stack is empty
			rebuild(); // thus we divide the left one into two
		}
		// now both stacks are non-empty
		if(stacks[1].top().second == start_pos + S.size() - k){
			stacks[1].pop();
		}		
	}else if(S.size() == k){ // stacks contain 1 element in total
		clear_stacks();
	}
	rem_right_update_S();
}	



void TwoStackMinimizer::add_left(char c){
	add_left_update_S(c);
	if(S.size() >= k){
		if((stacks[0].size() == 0) || (stacks[0].top().first >= lefthash)){
			stacks[0].push(make_pair(lefthash, start_pos));
		}
	}
}

void TwoStackMinimizer::rem_left() {
	if(S.size() > k){
		if(stacks[0].size() == 0){ // the left stack is empty
			rebuild(); // thus we divide the right one into two
		}
		// now both stacks are non-empty
		if(stacks[0].top().second == start_pos){
			stacks[0].pop();
		}		
	}else if(S.size() == k){ // stacks contain 1 element in total
		clear_stacks();
	}
	
	rem_left_update_S();
}	


INT TwoStackMinimizer::min_pos(){
	if(stacks[0].size() == 0){
		return stacks[1].top().second;
	}else if(stacks[1].size() == 0){
		return stacks[0].top().second;
	}else if(stacks[0].top().first <= stacks[1].top().first){
		return stacks[0].top().second;
	}else{
		return stacks[1].top().second;
	}
}

