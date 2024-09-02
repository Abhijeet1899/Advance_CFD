preconditioner_SIP: preconditioner_SIP.o residual.o matrix_vector_product.o vectorDotproduct.o SIP_Pinv_r.o two_norm.o
	gcc preconditioner_SIP.o residual.o matrix_vector_product.o vectorDotproduct.o SIP_Pinv_r.o two_norm.o -o preconditioner_SIP -lm
	
preconditioner_SIP.o: preconditioner_SIP.c head.h
	gcc -c preconditioner_ILU.c
	
SIP_Pinv_r.o: SIP_Pinv_r.c
	gcc -c SIP_Pinv_r.c
