#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

#include "MAT.h"
//字母说明，dec为分解，signal为信号，Len为长度，filter为滤波器，
 //EX为延拓，CON为卷积，DSam为下采样，D为下，sam为采样。L为低频信号，H为高频信号
 //ref为重构，
class TF
{
    public:
        TF(int dbn,int dbn_n,int row,int dataLen,int m,int n,int mode);
        ~TF();
        void TransForm_receive(pRWD pRWdate);                //得到读取的数据
        void TransForm_Matrix(pMAT pMatrix);                 //得到系数矩阵的数据
        void TransForm_print(int n);                         //打印分解函数
        void TransForm_row(int n);                           //获取单一通道数据
        void TransForm_EX(int n);                            //延拓函数
        void TransForm_CON(int n);                           //卷积函数
        void TransForm_DSam(int n);                          //下采样函数
        void transform(pRWD pRWdate,pMAT pMatrix);           //分解函数
        pDEC getTransform_DEC(){return pDEC_t;}              //返回分解结果 
        int getTransForm_DEC_Len(){return DEC_Len;}          //返回分解系数列数
        void TransForm_mkdir_file();                         //创建文件夹函数
        void TransForm_mkdir_txt(int n);                     //创建文件 

    private:   
        int DBN;                                             //DBn小波系数      
        int DBN_N;                                           //DBn小波分解层数
        int data_row;                                        //信号矩阵行长度
        int data_line;                                       //信号矩阵列长度
        int file_m;                                          //文件夹号
        int file_n;                                          //txt编号
        int mode;                                            //读取模式
        int dec_dataLen;                                     //分解信号长度
        int filterLen;                                       //滤波器长度
        int dec_EXLen;                                       //分解对称延拓长度
        int dec_CONLen;                                      //分解卷积长度
        int dec_DSamLen;                                     //分解下采样长度
        int DEC_Len;
        pRWD dec_pRWdate;                                    //接收数据
        pMAT pMatrix_t;                                      //接收矩阵数据
        pDEC pDEC_t;                                         //初始化分解矩阵
        std::string DIR;                                     //字符串
        std::vector<std::vector<double> > signal_row;        //一个行向量
        std::vector<std::vector<double> > dec_signal;        //分解信号
        std::vector<std::vector<double> > dec_EX;            //分解延拓信号    
        std::vector<std::vector<double> > dec_CON;           //分解卷积信号 
           
};


#endif