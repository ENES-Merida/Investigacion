#ifndef CONSTANTES_H
#define CONSTANTES_H

#include "herramientas.h"

const int mi = 6, nj = 5, nn = 1024;
const double pi = 3.1415926535;

#pragma acc routine vector
void ensambla_tdmax(
    double **AI,
    double **AC,
    double **AD,
    double **resultx,
    const double deltax,
    const double deltay,
    double **temp_ant,
    const double cond_ter,
    const double temp_ini,
    const double temp_fin,
    const int jj)
{
    int iin;

    /*
     * Se definen las condiciones de frontera
     */
    AC[0][jj] = 1.0;
    AD[0][jj] = 0.0;
    resultx[0][jj] = temp_ini;
    AI[mi - 1][jj] = 0.0;
    AC[mi - 1][jj] = 1.0;
    resultx[mi - 1][jj] = temp_fin;
/*
 * Ensamblado de la matriz tridiagonal y del vector de resultados
 */
#pragma acc loop vector
    for (iin = 1; iin < mi - 1; iin++)
    {
        AI[iin][jj] = -1.0 * cond_ter / (deltax * deltax);
        AC[iin][jj] = 2.0 * cond_ter * (1.0 / (deltax * deltax) + 1.0 / (deltay * deltay));
        AD[iin][jj] = -1.0 * cond_ter / (deltax * deltax);
        resultx[iin][jj] = cond_ter / (deltay * deltay) * temp_ant[iin][jj + 1] + cond_ter / (deltay * deltay) * temp_ant[iin][jj - 1];
    }
}

#pragma acc routine vector
void ensambla_tdmay(
    double **BI,
    double **BC,
    double **BD,
    double **resulty,
    const double deltax,
    const double deltay,
    double **temper,
    const double cond_ter,
    const double flux_aba,
    const double flux_arr,
    const int jj)
{
    int jjn;

    /*
     * Se definen las condiciones de frontera
     */
    BC[0][jj] = 1.0;             // -1.0 / deltay;
    BD[0][jj] = 0.0;             // 1.0 / deltay;
    resulty[0][jj] = 308.0;      // flux_aba;
    BI[nj - 1][jj] = 0.0;        // -1.0 / deltay;
    BC[nj - 1][jj] = 1.0;        // 1.0 / deltay;
    resulty[nj - 1][jj] = 308.0; // flux_arr
/*
 * Ensamblado de la matriz tridiagonal y del vector de resultados
 */
#pragma acc loop vector
    for (jjn = 1; jjn < nj - 1; jjn++)
    {
        BI[jjn][jj] = -1.0 * cond_ter / (deltay * deltay);
        BC[jjn][jj] = 2.0 * cond_ter * (1.0 / (deltay * deltay) + 1.0 / (deltax * deltax));
        BD[jjn][jj] = -1.0 * cond_ter / (deltay * deltay);
        resulty[jjn][jj] = cond_ter / (deltax * deltax) * temper[jj + 1][jjn] + cond_ter / (deltax * deltax) * temper[jj - 1][jjn];
    }
}

#pragma acc routine
void tri(
    double **a,
    double **b,
    double **c,
    double **r,
    const int filas,
    const int columna)
{
    int i;

    // Eliminacion elementos bajo la matriz - Eliminacion Gaussiana
    for (i = 1; i < filas; i++)
    {
        r[i][columna] = r[i][columna] - (a[i][columna] / b[i - 1][columna]) * r[i - 1][columna];
        b[i][columna] = b[i][columna] - (a[i][columna] / b[i - 1][columna]) * c[i - 1][columna];
    }
    // Solucion para r - Substitucion hacia atras
    r[filas - 1][columna] = r[filas - 1][columna] / b[filas - 1][columna];
    for (i = filas - 2; i >= 0; i--)
    {
        r[i][columna] = (r[i][columna] - c[i][columna] * r[i + 1][columna]) / b[i][columna];
    }
}

int obtener_nnz_matriz()
{
    int non_zero_elements = 12 + (mi + nj - 8) * 8 + (mi - 4) * (nj - 4) * 5;
    return non_zero_elements;
}

int obtener_indice_columna(int ii, int jj)
{
    int indice_columna;
    indice_columna = (jj - 1) * (mi - 2) + ii - (mi - 1) - 1;
    return indice_columna;
}

void band_matrix(double **BI, double **AI, double **AC, double **AD, double **BD, int nnz_e)
{
    double *csrVal;
    int *csrPtr;
    int *csrColInd;

    int size_ptr = (mi - 2) * (nj - 2) + 1;

    csrVal = allocate_memory_vector(nnz_e);
    csrColInd = allocate_memory_vector_int(nnz_e);
    csrPtr = allocate_memory_vector_int(size_ptr);

    int kk = -1;
    int tt = 0;
    int counter = -1;
    csrPtr[tt] = 0;

    for (int jj = 1; jj < nj - 1; jj++)
    {
        // printf("paso por aqui k=%d\n", kk);
        for (int ii = 1; ii < mi - 1; ii++)
        {
            if (jj >= 2 && jj < nj - 1)
            {
                kk++;
                csrVal[kk] = BI[jj - 1][ii];
                csrColInd[kk] = obtener_indice_columna(ii + 1, jj);
                counter++;
            }
            if (ii >= 2 && ii < mi - 1)
            {
                kk++;
                csrVal[kk] = AI[ii - 1][jj];
                csrColInd[kk] = obtener_indice_columna(ii, jj + 1);
                counter++;
            }
            kk++;
            if (kk >= nnz_e)
                break;
            // printf("paso por aqui k=%d\n", kk);
            csrVal[kk] = AC[ii][jj];
            csrColInd[kk] = obtener_indice_columna(ii + 1, jj + 1);
            counter++;
            if (ii >= 1 && ii < mi - 2)
            {
                kk++;
                csrVal[kk] = AD[ii + 1][jj];
                csrColInd[kk] = obtener_indice_columna(ii + 2, jj + 1);
                counter++;
            }
            if (jj >= 1 && jj < nj - 2)
            {
                kk++;
                csrVal[kk] = BD[jj + 1][ii];
                csrColInd[kk] = obtener_indice_columna(ii + 1, jj + 2);
                counter++;
            }
            tt++;
            csrPtr[tt] = counter + 1;
        }
    }
    // print_vector_int(csrColInd, nnz_e);
    // print_vector_int(csrPtr, (mi - 2) * (nj - 2) + 1);
    // printVector(csrVal, nnz_e);
    for (int ii = 0; ii < nnz_e; ii++)
    {
        if (ii < size_ptr)
        {
            printf("Val[%d] = %f | ColInd[%d] = %d | Ptr[%d] = %d\n", ii, csrVal[ii], ii, csrColInd[ii], ii, csrPtr[ii]);
            continue;
        }
        printf("Val[%d] = %f | ColInd[%d] = %d\n", ii, csrVal[ii], ii, csrColInd[ii]);
    }

    free(csrPtr);
    free(csrColInd);
    free(csrVal);
}

#endif