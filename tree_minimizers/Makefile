MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -fopenmp -msse4.2 -fomit-frame-pointer -funroll-loops 
 
LFLAGS= -std=c++11 -O3  -I .. -I ../common_code -Wno-deprecated
EXE=   tree_minimizers
 
SRC=    main.cpp ../common_code/input.cpp ../minimizer_queues.cpp minimizers.cpp minimizer_index.cpp ../common_code/krfp.cpp ../common_code/utils.cpp

OBJ2=    main.o input.o minimizer_queues.o minimizers.o minimizer_index.o krfp.o utils.o

HD=     ../common_code/input.h ../minimizer_queues.h minimizers.h minimizer_index.h ../common_code/krfp.h ../common_code/utils.h Makefile

# 
# No need to edit below this line 
# 
 
.SUFFIXES: 
.SUFFIXES: .cpp .o 
 
OBJ=    $(SRC:.cpp=.o) 
 
.cpp.o: 
	$(CC) $(CFLAGS)-c $(LFLAGS) $< 
 
all:    $(EXE) 
 
$(EXE): $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $(OBJ2) $(LFLAGS) 
 
$(OBJ): $(MF) $(HD) 
 
clean: 
	rm -f $(OBJ) $(OBJ2) $(EXE) *~

clean-all:
	rm -f $(OBJ) $(OBJ2) $(EXE) *~
	rm -r libsdsl
	rm -r sdsl-lite


