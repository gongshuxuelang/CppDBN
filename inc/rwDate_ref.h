#ifndef __RWDATE_REF_H__
#define __RWDATE_REF_H__

#include "MAT.h"

class RWDATE_REF
{
    public:
        RWDATE_REF(int dbn_n,int m,int n, int ft,pDEC pDEC_date);
        ~RWDATE_REF()
        {

        }
        bool ReadDate_ref();                            //读数据
        bool Write_Date_ref();                          //写数据   
        void Print_rwDate_ref();                        //打印数据 
        pREF getRaw_Data_ref(){return pRWdate_ref;}     //获得原始指针  
        int getFile_m_ref(){return file_m;}             //返回文件夹号
        int getFile_n_ref(){return file_n;}             //返回txt号  
        int getFile_txt_ref(){return file_txt;}         //返回txt文件      
        int getFile_data_line(){return data_line;}      //信号长度
    private:
        int file_m;                                     //读文件夹号
        int file_n;                                     //读文件夹中的txt
        int file_txt;                                   //行号
        int data_line;                                  //信号长度        
        int DBN_N;                                      //DBn小波系数        
        pREF pRWdate_ref;                               //原始数据
        pDEC pDEC_t;                                    //接收的分解后的结果
};

#endif