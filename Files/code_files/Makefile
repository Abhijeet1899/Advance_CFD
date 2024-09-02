#-*- MakeFile -*-

#target: dependencies
#	action

all: main preconditioner_jacobi preconditioner_ILU preconditioner_SIP

main: main.o residual.o matrix_vector_product.o vectorDotproduct.o two_norm.o
	gcc main.o residual.o matrix_vector_product.o vectorDotproduct.o two_norm.o -o main -lm

main.o: main.c head.h
	gcc -c main.c
	
matrix_vector_product.o: matrix_vector_product.c
	gcc -c matrix_vector_product.c
		
residual.o: residual.c
	gcc -c residual.c

vectorDotproduct.o: vectorDotproduct.c
	gcc -c vectorDotproduct.c
	
two_norm.o: two_norm.c
	gcc -c two_norm.c
	
preconditioner_jacobi: preconditioner_jacobi.o residual.o matrix_vector_product.o vectorDotproduct.o two_norm.o
	gcc preconditioner_jacobi.o residual.o matrix_vector_product.o vectorDotproduct.o two_norm.o -o preconditioner_jacobi -lm

preconditioner_jacobi.o: preconditioner_jacobi.c head.h
	gcc -c preconditioner_jacobi.c
	
preconditioner_ILU: preconditioner_ILU.o residual.o matrix_vector_product.o vectorDotproduct.o ILU_Pinv_r.o two_norm.o
	gcc preconditioner_ILU.o residual.o matrix_vector_product.o vectorDotproduct.o ILU_Pinv_r.o two_norm.o -o preconditioner_ILU -lm
	
preconditioner_ILU.o: preconditioner_ILU.c head.h
	gcc -c preconditioner_ILU.c
	
ILU_Pinv_r.o: ILU_Pinv_r.c
	gcc -c ILU_Pinv_r.c
	
preconditioner_SIP: preconditioner_SIP.o residual.o matrix_vector_product.o vectorDotproduct.o SIP_Pinv_r.o two_norm.o
	gcc preconditioner_SIP.o residual.o matrix_vector_product.o vectorDotproduct.o SIP_Pinv_r.o two_norm.o -o preconditioner_SIP -lm
	
preconditioner_SIP.o: preconditioner_SIP.c head.h
	gcc -c preconditioner_SIP.c
	
SIP_Pinv_r.o: SIP_Pinv_r.c
	gcc -c SIP_Pinv_r.c
	
