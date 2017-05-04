FLAGS = -lblas -lm -pg -O3 -std=gnu99 -L/usr/local/cuda/lib64 -lcudart

all: make run

make:
	nvcc -c dgemm-dlaiton.cu
	gcc main.c dgemm-dlaiton.o $(FLAGS) -o benchmark-dlaiton

run:
	./benchmark-dlaiton