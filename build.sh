mkdir -p bin
gcc src/test.c src/lib/sph.c src/lib/ker.c -o bin/test -lm -std=c2x
