#ifndef __REVERSE_TRANSFORM_H__
#define __REVERSE_TRANSFORM_H__

#include "MAT.h"
//字母说明，dec为分解，signal为信号，Len为长度，filter为滤波器，
//EX为延拓，CON为卷积，DSam为下采样，D为下，sam为采样。
//ref为重构，
//Alpha Beta  Delta Theta 
class RTF
{
    public:
        RTF(int dbn,int dbn_n,int dec_len,int fm,int fn,double HZ);
        ~RTF();
        void Reverse_transform_choose_signal(std::string str);            //波形选择
        void Reverse_transform_init(pMAT pMatrix,std::string str);        //重构初始化
        void Reverse_transform_Matrix(pMAT pMatrix);                      //接受系数矩阵
        void Reverse_transform_REF_DATA(pREF pREF_data);                  //接收分解的数据
        void Reverse_transform_choose();                                  //选择重构的波形
        void Reverse_transform_USam();                                    //重构上采样函数
        void Reverse_transform_EX();                                      //重构延拓函数
        void Reverse_transform_CON(int n);                                //重构卷积函数                  
        void Reverse_transform_print(int n);                              //打印重构结果
        void Reverse_transform(pREF pREF_data,int ft);                    //重构执行函数
        void Reverse_transForm_mkdir_file();                              //创建文件夹
        void Reverse_transForm_mkdir_txt(int n);                          //创建文件
    private:
        const double Alpha_begin = 7.81;                                  //Alpha起始频率
        const double Alpha_end   = 13.28;                                 //Alpha终止频率
        const double Beta_begin  = 13.28;                                 //Beta起始频率
        const double Beta_end    = 30.47;                                 //Beta终止频率
        const double Delta_begin = 0.78;                                  //Delta起始频率
        const double Delta_end   = 3.91;                                  //Delta终止频率    
        const double Theta_begin = 3.91;                                  //Theta起始频率    
        const double Theta_end   = 7.81;                                  //Theta终止频率
        int Alpha_row;                                                    //Alpha波的长度
        int Beta_row;                                                     //Beta波的长度
        int Delta_row;                                                    //Delta波的长度
        int Theta_row;                                                    //Theta波的长度
        int DBN;                                                          //分解DBN分解系数
        int DBN_N;                                                        //小波分解层数
        int REF_Len;                                                      //原始重构小波长度
        double frequency;                                                 //原始采样频率
        double tf;                                                        //频率步长
        int choose_row;
        double choose_begin;
        int file_m;
        int file_n;        
        int ref_dataLen;                                                  //重构小波长度      
        int filterLen;                                                    //滤波器长度 
        int ref_EXLen;                                                    //重构对称延拓长度
        int ref_CONLen;                                                   //重构卷积长度
        int ref_USamLen;                                                  //重构上采样长度
        
        std::string ref_str_signal;                                       //重构信号选择
        std::string DIR;
        pMAT pMatrix_t;                                                   //接收系数矩阵        
        pREF pREF_data_t;                                                 //存储重构的结果       
        std::vector<double> ref_signal;                                   //重构信号
        std::vector<double> ref_USam;                                     //重构上采样信号
        std::vector<double> ref_EX;                                       //重构延拓信号    
        std::vector<double> ref_CON;                                      //重构卷积信号 
        std::map<int,std::vector<double> > ref_signal_map;                //重构信号map
        std::map<int,std::vector<double> > ref_USam_map;                  //重构上采样map
        std::map<int,std::vector<double> > ref_EX_map;                    //重构延拓信号map
        std::map<int,std::vector<double> > ref_CON_map;                   //重构卷积信号map
        
};

#endif