The main directory contains implementation of the heap and two-stack structures from the paper (files minimizer_queues.cpp and minimizer_queues.h),
plus a simple use of such a structure (file main.cpp)

The remaining directories contains the code generating the input of the two applications shown in the paper as well as test the structures on those inputs.

Running both the window and trie experiments in the setting from the paper can be done using file experiments.sh
(first argument is the file with the data-set, second argument is z).

An example use:

bash experiments.sh input_data/efm.in 32
