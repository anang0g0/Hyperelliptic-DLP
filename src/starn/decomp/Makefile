all:
	clang -Wall -g -pg -O3  -mavx2 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp pat.c 

ber:
	gcc -Wall -g -pg -O3 -mavx2 -mtune=native -march=native -ffast-math -funroll-loops -fopenmp ber.c

euc:
	gcc -Wall -g -pg -O3 -mavx2 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp euc.c

pat:
	gcc -Wall -g -pg -mavx2 -O3 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp pat.c


clean:
	rm -f a.out

