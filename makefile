all: problem1 

problem1: Homework8.c
	gcc -o problem1 Homework8.c  -lm

clean:
	rm problem1 