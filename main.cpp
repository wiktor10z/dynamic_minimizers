#include <iostream>
#include "minimizer_queues.h"


int main(){
	HeapMinimizer hm(5);
	
	hm.add_right('A');
	hm.add_right('A');
	hm.add_right('C');
	hm.add_right('T');
	hm.add_right('G');
	hm.add_left('G');
	std::cout << hm.min_pos()<<endl;
	hm.rem_right();
	hm.rem_left();
	
}
