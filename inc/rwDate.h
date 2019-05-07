#ifndef __RWDATE_H__
#define __RWDATE_H__

#include "MAT.h"

class rwDate
{
    public:
        rwDate(int dbn_n,int row,int line,int m,int n, int mo)
        {
            DBN_N = dbn_n;
            data_row = row;
            data_line = line;
            file_m = m;
            file_n = n;
            mode = mo;
        }
        ~rwDate()
        {
            pRWdate -> Raw_Data.erase( pRWdate -> Raw_Data.begin(), pRWdate -> Raw_Data.end());
            delete pRWdate;
            pRWdate = NULL;
        }
        bool ReadDate();                     //读数据
        bool Write_Date();                   //写数据   
        void Print_rwDate();                 //打印数据 
        pRWD getRaw_Data(){return pRWdate;}  //获得原始指针  

    private:
        int file_m;                          //读文件夹号
        int file_n;                          //读文件夹中的txt
        int mode;                            //读取的模式
        int DBN_N;                           //DBn小波系数
        int data_row;                        //信号行
        int data_line;                       //信号列
        pRWD pRWdate;                        //原始数据

};



#endif