all: make run

make:
	g++ -O0 Dlaiton++.cpp -std=c++11 -o Dlaiton++

run: make
	./Dlaiton++

clean:
	rm -rf Dlaiton++