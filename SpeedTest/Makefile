TESTS = 100
ARQUIVO1 = soma.txt
ARQUIVO2 = multiplicacao.txt
ARQUIVO3 = ponto_flutuante.txt
ARQUIVO4 = multiplica_por1.txt
ARQUIVO5 = outside.txt
ARQUIVO6 = reverse.txt
ARQUIVO7 = salvar.txt

all: make run

make:
	cc -O1 speedTest.c -Wall -std=gnu99 -o speedTest

run: make
	./speedTest

run-full: make
	for((i=0;i<$(TESTS);i++)) do ./speedTest; done

clean:
	rm -rf speedTest $(ARQUIVO1) $(ARQUIVO2) $(ARQUIVO3) $(ARQUIVO4) $(ARQUIVO5) $(ARQUIVO6) $(ARQUIVO7)

#make++:
#	g++ -O1 speedTest.cpp -o speedTest
#
#run++: make++
#	./speedTest
#
#make++11:
#	g++ -O1 speedTest.cpp -std=c++11 -o speedTest
#
#run++11: make++11
#	./speedTest


