#include "head.h"
/*
    本次函数为读文件函数。
    读取文件函数有三个参数mode，这个参数的意义为读取哪一个文件
    0为读取EEG数据，1为读取label数据。
    参数 n是读那个文件，m是读那个文件夹 n,m从1开始，例如 n = 1,m = 2;读取的文件是s2_1
    读取的文件格式是sm_n mode= 0的时候。
    mode = 1的时候为 SLm；
*/
bool rwDate::ReadDate()
{
    pRWdate = new RWD; //申请原始数据内存；
    pRWdate -> Raw_Data = std::vector<std::vector<double> >(data_row, std::vector<double>(data_line));//c初始化内存大小   
    
    std::string file_directory  = "../date/EEG/s";//文件目录用来区分读取EEG和label的，EEG为s，label为SL
    std::string file_directory1 = "../date/label/SL";
    std::string file_temporary;//  这个是一个文件夹变量
    std::string A;              //这个是循环的得出的字符串
    std::string B;              //这个是循环的得出的字符串
    if(mode < 0 || mode > 1)
    {
        std::cout << "读取文件出错！" << std::endl;
        return false;
    }
    if(mode == 0)
    {
        //读取EEG信号       
        A = boost::lexical_cast<std::string>(file_m);//读取文件夹
        file_directory.append(A);//file_temporary = s1
        B = boost::lexical_cast<std::string>(file_n);//得出文件名
        file_temporary.append(file_directory);
        file_temporary.append("/s");
        file_temporary.append(A);
        file_temporary.append("_");
        file_temporary.append(B);
        file_temporary.append(".txt");
        std::cout << "file_temporary = " << file_temporary << std::endl;
        
        std::ifstream ifstr_data(file_temporary);	//读文件  
        //测试文件
        //std::ifstream ifstr_data("/home/lcx/code/testdate/1.txt");
        for(std::vector<double>::size_type i = 0; i != pRWdate -> Raw_Data.size(); ++i)//  把文件读到内存中
        {
            for(std::vector<double>::size_type j = 0; j != pRWdate -> Raw_Data[i].size(); ++j)
            {            
                ifstr_data >> pRWdate -> Raw_Data[i][j];                            
            }            
        }
    }                   
    else
    {        
        A = boost::lexical_cast<std::string>(file_m);//读取文件夹
        file_directory1.append(A);//file_temporary = SL1       
        file_directory1.append("/SL");
        file_directory1.append(A);  
        file_directory1.append(".txt");          
        std::cout << "file_temporary1 = " << file_directory1 << std::endl;       
    }   
    return true;
}
bool rwDate::Write_Date()
{
    return true;
}

void rwDate::Print_rwDate()
{
    std::cout << "print_rwDate:"<<  pRWdate -> Raw_Data.size() << std::endl;
    for(std::vector<double>::size_type i = 0; i != pRWdate -> Raw_Data.size(); ++i)
    {
        for(std::vector<double>::size_type j = 0; j != pRWdate -> Raw_Data[i].size(); ++j)
        {
            std::cout << pRWdate -> Raw_Data[i][j] << " ";
        }
        std::cout<< std::endl;
        std::cout << (pRWdate -> Raw_Data).size()<< std::endl;
        std::cout << (pRWdate -> Raw_Data[0]).size()<< std::endl;
    }
}