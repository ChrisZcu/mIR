CC = gcc
CXX = g++


TARGET = program
OBJ = advanced.o main.o algorithm.o  util.o\
	./celltree/cellTree.o\
 ./rtree/collection.o  ./rtree/filemem.o  ./rtree/global.o \
 ./rtree/hypercube.o ./rtree/param.o  ./rtree/point.o  ./rtree/rentry.o  ./rtree/rnode.o  ./rtree/rtree.o \
 ./rtree/skyline.o  ./rtree/tgs.o ./rtree/virtualRNode.o baseline.o

$(TARGET): $(OBJ)
	$(CXX) -g -o $@  $^  ./celltree/liblpsolve55.a libqhullcpp.a libqhullstatic_r.a -ldl -lm -no-pie


%.o: %.c
	$(CC) -g -c $< -o $@ 

%.o: %.cpp
	$(CXX) -g -c $< -o $@ 

.PHONY: clean
	
clean: 
	-rm -f *.o
	-rm -f $(TARGET)
	-rm -f ./celltree/*.o
	-rm -f ./rtree/*.o