//Inclue header file==============================================
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
//================================================================

//Matrix structure================================================
typedef struct _matrix
{
    int row_index;
    int col_index;
    float value;
} Matrix_t;


//Vector structure===============================================
typedef struct _vector
{
    int row_index;
    float value;
} Vector_t;


//Sparce Matrix=================================================
// typedef struct _s_matrix
// {
//     S_Matrix_t* next_cell;
//     float value;
// }S_Matrix_t;


//Create Matrix structure========================================
Matrix_t* create_matrix(int Num_Row, int Num_Col)
{
    int Num_cell = Num_Row * Num_Col;
    
    Matrix_t *New_Matrix = malloc(Num_cell * sizeof(Matrix_t));

    int cell_i = 0;
    for (int Col_i = 0; Col_i < Num_Col; ++Col_i)
    {
        for (int Row_i = 0; Row_i < Num_Row; ++Row_i)
            {
                New_Matrix[cell_i].col_index = Col_i;
                New_Matrix[cell_i].row_index = Row_i;
                New_Matrix[cell_i].value = 0.0;
                cell_i = cell_i + 1;
            }
        }
    return New_Matrix;
}


//Create Vector structure=======================================
Vector_t* create_vector(int Num_Row)
{
    Vector_t *New_Vector = malloc(Num_Row * sizeof(Vector_t));
    for (int Row_i = 0; Row_i < Num_Row; ++Row_i)
    {
        New_Vector[Row_i].row_index = Row_i;
        New_Vector[Row_i].value = 0.0;
    }
    return New_Vector;
}

//Function to convert Row-Col for 1 dimension index============
int convert_index(int Row_index, int Col_index, int s)
{
    return (Row_index * s) + Col_index;
}


//Function to calculate squere error===========================
float check_err(Matrix_t *M_, Vector_t *b_, Vector_t *x_, int s)
{
    float err = 0.0;
    for (int i = 0; i < s; ++i)
    {
        float Mx = 0;
        for (int j = 0; j < s; ++j)
        {
            Mx = Mx + M_[convert_index(i, j, s)].value * x_[j].value;
        }
        err = err + (b_[i].value - Mx)*(b_[i].value - Mx);
    }
    return err;
}


//Gausse sidel method==========================================
Vector_t* GS_method(Matrix_t* M_, Vector_t* b_, Vector_t* x_, int s)
{
    float err = 1.0;
    while(err>0.00001)
    {
        for(int i = 0; i<s; ++i)
        {
            float sum = 0.0;
            for(int j = 0; j<s; ++j)
            {
                if(j==i)
                {
                    1==1;
                }else{
                    sum = sum +\
                    M_[convert_index(i, j, s)].value * x_[j].value;
                }
            }
            x_[i].value = (1/M_[convert_index(i,i,s)].value)\
            *(b_[i].value - sum);
        }
        err = check_err(M_, b_, x_, s);
    }
    return x_;
}

//Function to check diagonaly dominant==========================
int check_diagonaly_dominant(Matrix_t* M_, Vector_t* b_, int s)
{
    for(int i = 0; i<s; ++i)
    {
        float sum = 0.0;
        for(int j = 0; j<s; ++j)
        {
            if(j==i)
            {
                1==1;
            }else{
                sum = sum +\
                abs(M_[convert_index(i, j, s)].value);
            }
        }
        if (abs(M_[convert_index(i, i, s)].value) <= sum)
        {
            return 1; 
        }
    }
    return 0;
}

//Check sparse matrix==========================
// S_Matrix_t cheak_sparce(Matrix_t *M_, int s)
// {
//     S_Matrix_t *SM = malloc(s * sizeof(S_Matrix_t));

//     for (int i = 0; i < s; ++i)
//     {
//         for(int j=0; j<s; ++j)
//         {
//             if(M_[convert_index(i,j,s)].value != 0)
//             {
//                 S_Matrix_t* NEW_SM = malloc(sizeof(S_Matrix_t));
                
//             }
//         }
//     }
// }


//Main function================================================
int main()
{
    int s = 3;
    Matrix_t* M = create_matrix(s, s);
    Vector_t* b = create_vector(s);
    Vector_t *x = create_vector(s);

    //Setting value
    M[convert_index(0, 0, s)].value = 3.0;
    M[convert_index(0, 1, s)].value = 1.0;
    M[convert_index(0, 2, s)].value = 1.0;
    b[0].value = 5.0;

    M[convert_index(1, 0, s)].value = 1.0;
    M[convert_index(1, 1, s)].value = -3.0;
    M[convert_index(1, 2, s)].value = 1.0;
    b[1].value = -1.0;

    M[convert_index(2, 0, s)].value = 2.0;
    M[convert_index(2, 1, s)].value = 2.0;
    M[convert_index(2, 2, s)].value = 5.0;
    b[2].value = 9.0;

    //Setting initial value
    x[0].value = 3.0;
    x[1].value = 3.0;
    x[2].value = 3.0;

    int flag = check_diagonaly_dominant(M,b,s);
    if(flag == 1)
    {
        printf("Not diagonally dominant matrix !\n");
        return 1;
    }

    
    //Gausse-Sidel
    GS_method(M, b, x, s);

    //Output solution
    for (int i = 0; i < s; ++i)
    {
        printf("x%i = %f\n", i, x[i].value);
    }
}
