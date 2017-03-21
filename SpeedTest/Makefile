all: make run

make:
	gcc -O0 speedTest.cpp -o speedTest

run: make
	./speedTest

clean:
	rm -rf speedTest

make++:
	g++ -O0 speedTest.cpp -o speedTest

run++: make++
	./speedTest

make++11:
	g++ -O0 speedTest.cpp -std=c++11 -o speedTest

run++11: make++11
	./speedTest