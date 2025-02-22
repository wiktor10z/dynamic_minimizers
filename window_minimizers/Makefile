MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -fopenmp -msse4.2 -fomit-frame-pointer -funroll-loops 
 
LFLAGS= -std=c++14 -O3 -DNDEBUG -I .. -I ../common_code

EXE=   window_minimizers 
 
SRC=    main.cpp ../common_code/input.cpp estimation.cpp  property_string.cpp minimizer_index.cpp ../common_code/krfp.cpp ../common_code/utils.cpp minimizers.cpp ../minimizer_queues.cpp

OBJ2=    main.o input.o estimation.o property_string.o minimizer_index.o  krfp.o utils.o minimizers.o minimizer_queues.o

HD=     ../common_code/input.h estimation.h property_string.h minimizer_index.h ../common_code/krfp.h ../common_code/utils.h minimizers.h ../minimizer_queues.h Makefile

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



