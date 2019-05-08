#include "head.h"

RWDATE_REF::RWDATE_REF(int dbn_n,int m,int n, int ft,pDEC pDEC_date)
{
    DBN_N = dbn_n;      //分解层数   
    file_m = m;         //文件夹号
    file_n = n;         //子文件夹号
    file_txt = ft;     //txt号
    pDEC_t = pDEC_date;
    data_line = pDEC_t -> DEC_date[0].size();
}

bool RWDATE_REF::ReadDate_ref() //读数据
{
    pRWdate_ref = new REF; //申请原始数据内存；
    pRWdate_ref -> REF_date = std::vector<std::vector<double> >(power(2,DBN_N), std::vector<double>(data_line));//c初始化内存大小   
    
    std::string file_temporary  = "../TF_DATE/s";//文件目录用来区分读取EEG和label的，EEG为s，
   
    std::string A;              //这个是循环的得出的字符串
    std::string B;              //这个是循环的得出的字符串
    std::string C;              //这个是循环的得出的字符串
    //读取分解信号       
    A = boost::lexical_cast<std::string>(file_m);//读取文件夹
    B = boost::lexical_cast<std::string>(file_n);//得出文件名
    C = boost::lexical_cast<std::string>(file_txt);//得出文件名
    file_temporary = file_temporary + A + "/s" + A + "_" + B + "/s" + A + "_" + B + "_" + C +".txt";
    std::cout << "file_temporary = " << file_temporary << std::endl;
    
    std::ifstream ifstr_data(file_temporary);	//读文件  

    for(std::vector<double>::size_type i = 0; i != pRWdate_ref -> REF_date.size(); ++i)//  把文件读到内存中
    {
        for(std::vector<double>::size_type j = 0; j != pRWdate_ref -> REF_date[i].size(); ++j)
        {            
            ifstr_data >> pRWdate_ref -> REF_date[i][j];                            
        }            
    }
    return true;
}

void RWDATE_REF::Print_rwDate_ref()//打印数据 
{
    std::cout << "重构文件" << std::endl;
    for(std::vector<double>::size_type i = 0; i != pRWdate_ref -> REF_date.size(); ++i)
    {
        for(std::vector<double>::size_type j = 0; j != pRWdate_ref -> REF_date[i].size(); ++j)
        {
            std::cout << pRWdate_ref -> REF_date[i][j] << " ";
        }
    }
}