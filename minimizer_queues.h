#include <iostream>
#include <deque>
#include <set>
#include <stack>
#include "krfp.h"

using namespace std;

typedef uint64_t UINT;
typedef int64_t INT;

class MinimizerQueue {
protected:
	INT start_pos = 0;
	UINT k,lefthash = 0, righthash = 0;
	deque<char> S;
	
	void add_right_update_S(char c);
	void rem_right_update_S();
	void add_left_update_S(char c);
	void rem_left_update_S();
public:
	MinimizerQueue(UINT k1) : k(k1){
		karp_rabin_hashing::init(k); // initialization with fixed k
	}
	virtual void add_right(char c) = 0;
	virtual void rem_right() = 0;
	virtual void add_left(char c) = 0;
	virtual void rem_left() = 0;
	virtual INT min_pos() = 0;
	
	void move_left(char c){
		add_left(c);
		rem_right();
	}
	
	void move_right(char c){
		add_right(c);
		rem_left();
	}	
};


class HeapMinimizer : public MinimizerQueue {
private:
	set<pair<UINT, INT>> heap; // (KR hash, position)
public:
	HeapMinimizer(UINT k1) : MinimizerQueue(k1){}
	void add_right(char c) override;
	void rem_right() override;
	void add_left(char c) override;
	void rem_left() override;
	INT min_pos() override;
};


class TwoStackMinimizer : public MinimizerQueue {
private:
	stack<pair<UINT, INT>> stacks[2]; // (KR hash, position) first one stores left part, second one stores the right part
	void clear_stacks();
	void rebuild();

public:
	TwoStackMinimizer(UINT k1) : MinimizerQueue(k1){}
	void add_right(char c) override;
	void rem_right() override;
	void add_left(char c) override;
	void rem_left() override;
	INT min_pos() override;
};

