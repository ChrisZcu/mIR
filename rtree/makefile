CC = gcc
CXX = g++


TARGET = program
OBJ = virtualRNode.o collection.o filemem.o global.o hypercube.o main.o param.o point.o rentry.o rnode.o rtree.o skyline.o tgs.o
	
$(TARGET): $(OBJ)
	$(CXX) -o $@ $^


%.o: %.c
	$(CC)  -c $< -o $@ 

%.o: %.cpp
	$(CXX) -c $< -o $@

.PHONY: clean
	
clean: 
	-rm -f *.o
	-rm -f $(TARGET)