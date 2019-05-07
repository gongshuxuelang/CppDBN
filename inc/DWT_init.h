#ifndef __DWT_INIT_H__
#define __DWT_INIT_H__
/*
此类用于初始化小波分解所有初始化的数据
比如滤波器的长度，信号的长度，信号的宽度，

*/
class DWT
{
    public:
        DWT(int dbn,int dbn_n,int row,int line,int mo,int fm,int fn,int ft)
        {
            DBN       = dbn;
            DBN_N     = dbn_n;
            data_row  = row;
            data_line = line;
            file_m    = fm;
            file_n    = fn;
            mode      = mo;
            file_txt  = ft;
           
        }    
        ~DWT(){}
        int getDWT_DBN(){return DBN;}               //获取DBN_N 值
        int getDWT_DBN_N(){return DBN_N;}           //获取DBN_N 值
        int getDWT_data_row(){return data_row;}     //获取data_row值
        int getDWT_data_line(){return data_line;}   //获取data_line值
        int getDWT_file_m(){return file_m;}         //获取file_m
        int getDWT_file_n(){return file_n;}         //获取file_n
        int getDWT_mode(){return mode;}             //获取模式
        int getDWT_file_txt(){return file_txt;}
        
    private:   
        int DBN;                                    //小波分解DBn系数
        int DBN_N;                                  //DBn小波层数
        int file_m;                                 //m的文件夹
        int file_n;                                 //每个文件夹有n个文件
        int mode;                                   //模式
        int file_txt;
        int data_row;                               //信号矩阵行长度
        int data_line;                              //信号矩阵列长度
               
};


#endif