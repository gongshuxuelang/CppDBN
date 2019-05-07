#include "head.h"

void Matrix::matrix()
{
    Matrix_Init(); //初始化系数矩阵
    creatMatrix();  //创建矩阵
}
bool Matrix::Matrix_Init()
{	
	pMatrix = new MAT;  
	pMatrix -> filter = new double*[4];
    for (int i = 0; i < 4; ++i)
	{		
        pMatrix -> filter[i] = new double[2 * DBN_N];
    }
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 2 * DBN_N; ++j)
		{
			pMatrix->filter[i][j] = 0;
		}
	}
	return true;
}

bool Matrix::creatMatrix()
{
	if (0 > DBN_N || DBN_N > 10)
	{
		std::cout << "你输入的维度有误！" << std::endl;	
		exit(-1);
	}
	
    switch (DBN_N)
	{
        case 1:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix -> filter[0][i] = db1_Lo_D[2 * DBN_N - 1 - i];
                pMatrix -> filter[1][i] = db1_Hi_D[2 * DBN_N - 1 - i];
                pMatrix -> filter[2][i] = db1_Lo_R[2 * DBN_N - 1 - i];
                pMatrix -> filter[3][i] = db1_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 2:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db2_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db2_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db2_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db2_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 3:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db3_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db3_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db3_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db3_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 4:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db4_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db4_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db4_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db4_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 5:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db5_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db5_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db5_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db5_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 6:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db6_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db6_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db6_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db6_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 7:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db7_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db7_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db7_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db7_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 8:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db8_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db8_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db8_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db8_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 9:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db9_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db9_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db9_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db9_Hi_R[2 * DBN_N - 1 - i];
            }
            break;
        case 10:
            for (int i = 0; i < 2 * DBN_N; ++i)
            {
                pMatrix->filter[0][i] = db10_Lo_D[2 * DBN_N - 1 - i];
                pMatrix->filter[1][i] = db10_Hi_D[2 * DBN_N - 1 - i];
                pMatrix->filter[2][i] = db10_Lo_R[2 * DBN_N - 1 - i];
                pMatrix->filter[3][i] = db10_Hi_R[2 * DBN_N - 1 - i];
            }
		    break;
	}
    return true;
}

 void Matrix::Print_matrix()
 {
     for(int i = 0; i < 4; ++i)
     {
         for(int j = 0; j < 2 * DBN_N; ++j)
         {
             std::cout << pMatrix -> filter[i][j] << "    ";
         }
         std::cout << std::endl;
     }
 }