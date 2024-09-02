#-*- MakeFile -*-

#target: dependencies
#	action

#all: 1 2

v: v.o v_cycle.o residual.o ristriction.o gauss_iteration.o prolongation.o two_norm.o
	gcc v.o v_cycle.o residual.o ristriction.o gauss_iteration.o prolongation.o two_norm.o -o v -lm

v.o: v.c head.h
	gcc -c v.c
	
v_cycle.o: v_cycle.c
	gcc -c v_cycle.c
		
residual.o: residual.c
	gcc -c residual.c

ristriction.o: ristriction.c
	gcc -c ristriction.c

guass_iteration.o: gauss_iteration.c
	gcc -c gauss_iteration.c

prolongation.o: prolongation.c
	gcc -c prolongation.c
	
2_norm.o: two_norm.c
	gcc -c two_norm.c
