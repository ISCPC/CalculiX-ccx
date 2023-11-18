#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include "matrix_io.h"

int save_matrix_csr(char *filename, INT_T neq, INT_T ndim, INT_T* pointers, INT_T* indice, double* value, INT_T flags) {
    int fd;
    INT_T info[8];
    size_t datasz;

    if ((fd = open(filename, O_CREAT|O_RDWR|O_TRUNC, 0644)) < 0) {
        return -1;
    }

    printf("INFO: Writing Matrix-A neq=%d, ndim=%d\n", neq, ndim);

    /* 
     * Put data
     */
    info[0] = neq;
    info[1] = ndim;
    info[2] = 0;
    //info[3] = SPMATRIX_TYPE_CSC | SPMATRIX_TYPE_INDEX0 | SPMATRIX_TYPE_SYMMETRIC;
    info[3] = flags;
    info[4] = neq;

    datasz = sizeof(INT_T)*8;
    if (write(fd, info, datasz) != datasz) {
        return -1;
    }

    datasz = sizeof(INT_T)*(neq+1);
    if (write(fd, pointers, datasz) != datasz) {
        printf("ERROR: cannot save pointers of A from %s\n", filename);
        return -1;
    }

    datasz = sizeof(INT_T)*ndim;
    if (write(fd, indice, datasz) != datasz) {
        printf("ERROR: cannot save indice of A from %s\n", filename);
        return -1;
    }

    datasz = sizeof(double)*ndim;
    if (write(fd, value, datasz) != datasz) {
        printf("ERROR: cannot save values of A from %s\n", filename);
        return -1;
    }

    close(fd);

    return 0;
}

int save_vector(char *filename, INT_T size, double* v) {
    int fd;
    size_t datasz;

    printf("INFO: Writing Vector size=%d\n", size);

    if ((fd = open(filename, O_CREAT|O_RDWR|O_TRUNC, 0644)) < 0) {
        return -1;
    }

    datasz = sizeof(INT_T);
    if (write(fd, &size, datasz) != datasz) {
        return -1;
    }

    datasz = sizeof(double)*size;
    if (write(fd, v, datasz) != datasz) {
        printf("ERROR: cannot save vector from %s\n", filename);
        return -1;
    }

    close(fd);

    return 0;
}



