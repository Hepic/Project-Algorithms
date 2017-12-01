CC = g++
FLAGS = -g -Wall
OUT = run

OBJS = main.o file_functions.o curve.o help_functions.o cluster.o distances.o  binary_mean_tree.o

run: $(OBJS)
	$(CC) $(FLAGS) $^ -o $(OUT)

main.o: main.cpp
	$(CC) $(FLAGS) -c $< -o $@

%.o: %.cpp %.h
	$(CC) $(FLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(OUT) 
