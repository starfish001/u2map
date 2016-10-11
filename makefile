CC	= gcc  -Wall

u2map:u2map.o group.o order.o prepare.o correct.o calculate.o matrix.o
	$(CC) -lm -o u2map u2map.o group.o order.o prepare.o correct.o calculate.o matrix.o
u2map.o:u2map.c prepare.h order.h group.h
	$(CC) -c u2map.c
group.o:group.c group.h order.h
	$(CC) -c group.c
order.o:order.c order.h
	$(CC) -c order.c
prepare.o:prepare.c prepare.h correct.h calculate.h
	$(CC) -c prepare.c
correct.o:correct.c correct.h calculate.h 
	$(CC) -c correct.c
calculate.o:calculate.c calculate.h matrix.h
	$(CC) -c calculate.c
matrix.o:matrix.c matrix.h
	$(CC) -c matrix.c
clean:
	rm u2map.o group.o order.o prepare.o correct.o calculate.o matrix.o
