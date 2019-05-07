#ifndef __MAT_H__
#define __MAT_H__



typedef struct MATRIX			//创建系数矩阵
{    
	double** filter;
}MAT, *pMAT;

typedef struct RWDATE// 二位数据组，读取数据矩阵
{    
	std::vector<std::vector<double> >  Raw_Data;
}RWD, *pRWD;

typedef struct DEC_DATA// 二位数据组，分解的矩阵
{ 
	std::vector<std::vector<double> >  DEC_date;
}DEC, *pDEC;

typedef struct REF_DATA// 二位数据组，重构的矩阵
{    
	std::vector<std::vector<double> >  REF_date;
}REF, *pREF;


#endif
