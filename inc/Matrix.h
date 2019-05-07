#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "MAT.h"

class Matrix
{
    public:
        Matrix(int dbn_n):DBN_N(dbn_n){}
        ~Matrix()
        {
           
        }
        bool Matrix_Init();                     //初始化系数矩阵
        bool creatMatrix();                     //创建系数矩阵
        void Print_matrix();                    //打印矩阵
        void matrix();                          //矩阵主函数               
        pMAT getMatrix(){return pMatrix;}       //获得系数矩阵
    private:
        int DBN_N;                              //DBN小波
        pMAT pMatrix;                           //系数矩阵
};

#endif