MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -fopenmp -msse4.2 -fomit-frame-pointer -funroll-loops 
 
LFLAGS= -std=c++20 -O3 -DNDEBUG -lboost_iostreams  -I common_code

EXE=   minimizer_queue
 
SRC=    main.cpp minimizer_queues.cpp common_code/krfp.cpp common_code/utils.cpp

OBJ2=   main.o minimizer_queues.o krfp.o utils.o

HD=     minimizer_queues.h common_code/krfp.h common_code/utils.h Makefile

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
