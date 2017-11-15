CC = g++
FLAGS = -g -c
OUT = run
OBJS = main.o file_functions.o curve.o help_functions.o

run: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT)

%.o: %.cpp %.h
	$(CC) $(FLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(OUT) 
