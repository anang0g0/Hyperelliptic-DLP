all:
	gcc -Wall -g -pg -O3  -mavx2 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp main.c

amd:
	gcc -Wall -g -pg -O3 -mavx2 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp main.c

ber:
	gcc -Wall -g -pg -O3 -mavx2 -mtune=native -march=native -ffast-math -funroll-loops -fopenmp Berlekamp.c

clang:
	clang -Wall -g -pg -O3 -mavx2 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp main.c

euc:
	gcc -Wall -g -pg -O3 -mavx2 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp Xeuclide.c

pat:
	gcc -Wall -g -pg -mavx2 -O3 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp Patterson.c


clean:
	rm -f a.out

