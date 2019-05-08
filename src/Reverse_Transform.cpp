#include "head.h"

RTF::RTF(int dbn,int dbn_n,int dec_len,int fm,int fn,double HZ):DBN(dbn),DBN_N(dbn_n),REF_Len(dec_len),frequency(HZ),file_m(fm),file_n(fn)
{
    filterLen   = 2 * DBN;                          //滤波器长度
    ref_USamLen = 2 * REF_Len + 1;                  //重构长采样长度
    ref_EXLen   = 2 * REF_Len + 2 * filterLen -1;   //重构对称延拓长度
    ref_CONLen  = 2 * REF_Len + filterLen;          //重构卷积长度
    ref_dataLen = 2 * REF_Len - filterLen + 2;      //重构信号长度
    tf = (frequency / 2) / (static_cast<double>(power(2,DBN_N)));//频率步长  
    Alpha_row = ceil(Alpha_end / tf) - ceil(Alpha_begin / tf) + 1;//计算重构信号量的行数 ceil()向上取整floor()向下取整
    Beta_row  = ceil(Beta_end  / tf) - ceil(Beta_begin  / tf) + 1;
    Delta_row = ceil(Delta_end / tf) - ceil(Delta_begin / tf) + 1;
    Theta_row = ceil(Theta_end / tf) - ceil(Theta_begin / tf) + 1;
}
RTF::~RTF()
{   

}
void RTF::Reverse_transform_choose_signal(std::string str)
{ 
    ref_str_signal = str;
    if("Alpha" == ref_str_signal)
    {
        choose_row   = Alpha_row;
        choose_begin = Alpha_begin;
    }else if("Beta" == ref_str_signal)
    {
        choose_row   = Beta_row;
        choose_begin = Beta_begin;
    }else if("Delta" == ref_str_signal)
    {   
        choose_row   = Delta_row;
        choose_begin = Delta_begin;
    }else if("Theta" == ref_str_signal)
    {
        choose_row   = Theta_row;
        choose_begin = Theta_begin;
    }
    return;
}
void RTF::Reverse_transform_init(pMAT pMatrix,std::string str)
{
    Reverse_transform_Matrix(pMatrix);          //接收矩阵    
    Reverse_transform_choose_signal(str);       //接收波形 
    
}
void RTF::Reverse_transForm_mkdir_file()
{
    DIR = "../RTF_DATE/" + ref_str_signal +"/s";
    std::string FILE_NAME = boost::lexical_cast<std::string>(file_m);
    DIR = DIR + FILE_NAME;
    DIR = "mkdir -p " + DIR;
    system(DIR.c_str());//创建Sn总文件夹
    std::string  num= boost::lexical_cast<std::string>(file_n);
    DIR = DIR + "/s"+ FILE_NAME+ "_" + num;
    system(DIR.c_str());//创建Sn_n总文件夹
    std::cout << DIR << std::endl;
}  
void RTF::Reverse_transForm_mkdir_txt(int n)
{
    DIR = "../RTF_DATE/"+ ref_str_signal +"/s";
    std::string FILE_NAME = boost::lexical_cast<std::string>(file_m);
    std::string  num= boost::lexical_cast<std::string>(file_n);
    std::string  NEMBER_NAME= boost::lexical_cast<std::string>(n);
    DIR = DIR + FILE_NAME;// /home/lcx/code/REF_DATE/S1
    DIR = DIR + "/s"+ FILE_NAME+ "_" + num + "/s" + FILE_NAME + "_" + num + "_" + NEMBER_NAME + ref_str_signal + ".txt"; //  /s1/s1_1/s1_1_1_Alpha.txt    
    std::cout << DIR << std::endl;
}
void RTF::Reverse_transform(pREF pREF_data,int ft)
{
    Reverse_transform_REF_DATA(pREF_data);      //接收重构信号
    Reverse_transform_choose();  
    Reverse_transForm_mkdir_file();
    Reverse_transForm_mkdir_txt(ft);               //根据波形选择数据

    for(int i = 0; i < DBN_N; ++i)
    {
        std::cout << "第"<< i << "重构" <<std::endl;
        Reverse_transform_USam();                   //数据上采样
        Reverse_transform_EX();
        Reverse_transform_CON(i);
        
        //重新计算系数
        ref_USamLen = 2 * REF_Len + 1;                  //重构长采样长度
        ref_EXLen   = 2 * REF_Len + 2 * filterLen -1;   //重构对称延拓长度
        ref_CONLen  = 2 * REF_Len + filterLen;          //重构卷积长度
        ref_dataLen = 2 * REF_Len - filterLen + 2;      //重构信号长度     
    }
    //Reverse_transform_print(4);
}

void RTF::Reverse_transform_choose()//重构波形选择,把vector封装成map
{
    ref_signal = std::vector<double>(REF_Len);      
    for(std::vector<double>::size_type i = 0; i < static_cast<std::vector<double>::size_type>(choose_row); ++i)
    {
        for(std::vector<double>::size_type j = 0; j < static_cast<std::vector<double>::size_type>(REF_Len); ++j)
        {
            ref_signal[j] = pREF_data_t -> REF_date[i + choose_begin / tf][j];                
        }        
        ref_signal_map.insert(make_pair(power(2,DBN_N) + choose_begin / tf + i, ref_signal));
    }  
    return;
}

void RTF::Reverse_transform_CON(int n)
{  
    ref_CON = std::vector<double>(ref_CONLen);          
    std::map<int,std::vector<double> >::iterator I = ref_EX_map.begin(); 
    std::map<int,std::vector<double> >::size_type Iend = ref_EX_map.size(); 

    int chn = floor(ref_EX_map.size() / 2);//计算上一层的个数的
    if( I -> first % 2 != 0 && (Iend + I -> first - 1) % 2 == 0) 
    {
        chn++;
    }   
    std::map<int,std::vector<double> >::iterator i_1 = I;//用于卷积的
    std::map<int,std::vector<double> >::iterator i_2 = (I++); 
    switch(DBN)
    {
        case 1:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *j_1 * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *j_2 * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1]; 
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *j_1 * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *j_1 * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1];
                    } 
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;    
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 2:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3]; 
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 3:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3] +
                                         *(j_2 + 4) * pMatrix_t -> filter[3][4] + *(j_2 + 5) * pMatrix_t -> filter[3][5];  
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3] +
                                        *(j_1 + 4) * pMatrix_t -> filter[3][4] + *(j_1 + 5) * pMatrix_t -> filter[3][5]; 
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 4:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3] +
                                         *(j_2 + 4) * pMatrix_t -> filter[3][4] + *(j_2 + 5) * pMatrix_t -> filter[3][5] +
                                         *(j_2 + 6) * pMatrix_t -> filter[3][6] + *(j_2 + 7) * pMatrix_t -> filter[3][7];  
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3] +
                                        *(j_1 + 4) * pMatrix_t -> filter[3][4] + *(j_1 + 5) * pMatrix_t -> filter[3][5] +
                                        *(j_1 + 6) * pMatrix_t -> filter[3][6] + *(j_1 + 7) * pMatrix_t -> filter[3][7];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 5:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3] +
                                         *(j_2 + 4) * pMatrix_t -> filter[3][4] + *(j_2 + 5) * pMatrix_t -> filter[3][5] +
                                         *(j_2 + 6) * pMatrix_t -> filter[3][6] + *(j_2 + 7) * pMatrix_t -> filter[3][7] +
                                         *(j_2 + 8) * pMatrix_t -> filter[3][8] + *(j_2 + 9) * pMatrix_t -> filter[3][9];   
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3] +
                                        *(j_1 + 4) * pMatrix_t -> filter[3][4] + *(j_1 + 5) * pMatrix_t -> filter[3][5] +
                                        *(j_1 + 6) * pMatrix_t -> filter[3][6] + *(j_1 + 7) * pMatrix_t -> filter[3][7] +
                                        *(j_1 + 8) * pMatrix_t -> filter[3][8] + *(j_1 + 9) * pMatrix_t -> filter[3][9];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 6:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3] +
                                         *(j_2 + 4) * pMatrix_t -> filter[3][4] + *(j_2 + 5) * pMatrix_t -> filter[3][5] +
                                         *(j_2 + 6) * pMatrix_t -> filter[3][6] + *(j_2 + 7) * pMatrix_t -> filter[3][7] +
                                         *(j_2 + 8) * pMatrix_t -> filter[3][8] + *(j_2 + 9) * pMatrix_t -> filter[3][9] +
                                         *(j_2 + 10) * pMatrix_t -> filter[3][10] + *(j_2 + 11) * pMatrix_t -> filter[3][11];     
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3] +
                                        *(j_1 + 4) * pMatrix_t -> filter[3][4] + *(j_1 + 5) * pMatrix_t -> filter[3][5] +
                                        *(j_1 + 6) * pMatrix_t -> filter[3][6] + *(j_1 + 7) * pMatrix_t -> filter[3][7] +
                                        *(j_1 + 8) * pMatrix_t -> filter[3][8] + *(j_1 + 9) * pMatrix_t -> filter[3][9] +
                                        *(j_1 + 10) * pMatrix_t -> filter[3][10] + *(j_1 + 11) * pMatrix_t -> filter[3][11];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 7:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_1 + 12) * pMatrix_t -> filter[2][12] + *(j_1 + 13) * pMatrix_t -> filter[2][13] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3] +
                                         *(j_2 + 4) * pMatrix_t -> filter[3][4] + *(j_2 + 5) * pMatrix_t -> filter[3][5] +
                                         *(j_2 + 6) * pMatrix_t -> filter[3][6] + *(j_2 + 7) * pMatrix_t -> filter[3][7] +
                                         *(j_2 + 8) * pMatrix_t -> filter[3][8] + *(j_2 + 9) * pMatrix_t -> filter[3][9] +
                                         *(j_2 + 10) * pMatrix_t -> filter[3][10] + *(j_2 + 11) * pMatrix_t -> filter[3][11] +
                                         *(j_2 + 12) * pMatrix_t -> filter[3][12] + *(j_2 + 13) * pMatrix_t -> filter[3][13];      
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3] +
                                        *(j_1 + 4) * pMatrix_t -> filter[3][4] + *(j_1 + 5) * pMatrix_t -> filter[3][5] +
                                        *(j_1 + 6) * pMatrix_t -> filter[3][6] + *(j_1 + 7) * pMatrix_t -> filter[3][7] +
                                        *(j_1 + 8) * pMatrix_t -> filter[3][8] + *(j_1 + 9) * pMatrix_t -> filter[3][9] +
                                        *(j_1 + 10) * pMatrix_t -> filter[3][10] + *(j_1 + 11) * pMatrix_t -> filter[3][11] +
                                        *(j_1 + 12) * pMatrix_t -> filter[3][12] + *(j_1 + 13) * pMatrix_t -> filter[3][13];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_1 + 12) * pMatrix_t -> filter[2][12] + *(j_1 + 13) * pMatrix_t -> filter[2][13];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 8:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_1 + 12) * pMatrix_t -> filter[2][12] + *(j_1 + 13) * pMatrix_t -> filter[2][13] +
                                         *(j_1 + 14) * pMatrix_t -> filter[2][14] + *(j_1 + 15) * pMatrix_t -> filter[2][15] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3] +
                                         *(j_2 + 4) * pMatrix_t -> filter[3][4] + *(j_2 + 5) * pMatrix_t -> filter[3][5] +
                                         *(j_2 + 6) * pMatrix_t -> filter[3][6] + *(j_2 + 7) * pMatrix_t -> filter[3][7] +
                                         *(j_2 + 8) * pMatrix_t -> filter[3][8] + *(j_2 + 9) * pMatrix_t -> filter[3][9] +
                                         *(j_2 + 10) * pMatrix_t -> filter[3][10] + *(j_2 + 11) * pMatrix_t -> filter[3][11] +
                                         *(j_2 + 12) * pMatrix_t -> filter[3][12] + *(j_2 + 13) * pMatrix_t -> filter[3][13] +
                                         *(j_2 + 14) * pMatrix_t -> filter[3][14] + *(j_2 + 15) * pMatrix_t -> filter[3][15];      
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3] +
                                        *(j_1 + 4) * pMatrix_t -> filter[3][4] + *(j_1 + 5) * pMatrix_t -> filter[3][5] +
                                        *(j_1 + 6) * pMatrix_t -> filter[3][6] + *(j_1 + 7) * pMatrix_t -> filter[3][7] +
                                        *(j_1 + 8) * pMatrix_t -> filter[3][8] + *(j_1 + 9) * pMatrix_t -> filter[3][9] +
                                        *(j_1 + 10) * pMatrix_t -> filter[3][10] + *(j_1 + 11) * pMatrix_t -> filter[3][11] +
                                        *(j_1 + 12) * pMatrix_t -> filter[3][12] + *(j_1 + 13) * pMatrix_t -> filter[3][13] +
                                        *(j_1 + 14) * pMatrix_t -> filter[3][14] + *(j_1 + 15) * pMatrix_t -> filter[3][15];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_1 + 12) * pMatrix_t -> filter[2][12] + *(j_1 + 13) * pMatrix_t -> filter[2][13] +
                                         *(j_1 + 14) * pMatrix_t -> filter[2][14] + *(j_1 + 15) * pMatrix_t -> filter[2][15];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 9:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_1 + 12) * pMatrix_t -> filter[2][12] + *(j_1 + 13) * pMatrix_t -> filter[2][13] +
                                         *(j_1 + 14) * pMatrix_t -> filter[2][14] + *(j_1 + 15) * pMatrix_t -> filter[2][15] +
                                         *(j_1 + 16) * pMatrix_t -> filter[2][16] + *(j_1 + 17) * pMatrix_t -> filter[2][17] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3] +
                                         *(j_2 + 4) * pMatrix_t -> filter[3][4] + *(j_2 + 5) * pMatrix_t -> filter[3][5] +
                                         *(j_2 + 6) * pMatrix_t -> filter[3][6] + *(j_2 + 7) * pMatrix_t -> filter[3][7] +
                                         *(j_2 + 8) * pMatrix_t -> filter[3][8] + *(j_2 + 9) * pMatrix_t -> filter[3][9] +
                                         *(j_2 + 10) * pMatrix_t -> filter[3][10] + *(j_2 + 11) * pMatrix_t -> filter[3][11] +
                                         *(j_2 + 12) * pMatrix_t -> filter[3][12] + *(j_2 + 13) * pMatrix_t -> filter[3][13] +
                                         *(j_2 + 14) * pMatrix_t -> filter[3][14] + *(j_2 + 15) * pMatrix_t -> filter[3][15] +
                                         *(j_2 + 16) * pMatrix_t -> filter[3][16] + *(j_2 + 17) * pMatrix_t -> filter[3][17];        
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3] +
                                        *(j_1 + 4) * pMatrix_t -> filter[3][4] + *(j_1 + 5) * pMatrix_t -> filter[3][5] +
                                        *(j_1 + 6) * pMatrix_t -> filter[3][6] + *(j_1 + 7) * pMatrix_t -> filter[3][7] +
                                        *(j_1 + 8) * pMatrix_t -> filter[3][8] + *(j_1 + 9) * pMatrix_t -> filter[3][9] +
                                        *(j_1 + 10) * pMatrix_t -> filter[3][10] + *(j_1 + 11) * pMatrix_t -> filter[3][11] +
                                        *(j_1 + 12) * pMatrix_t -> filter[3][12] + *(j_1 + 13) * pMatrix_t -> filter[3][13] +
                                        *(j_1 + 14) * pMatrix_t -> filter[3][14] + *(j_1 + 15) * pMatrix_t -> filter[3][15] +
                                        *(j_1 + 16) * pMatrix_t -> filter[3][16] + *(j_1 + 17) * pMatrix_t -> filter[3][17];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_1 + 12) * pMatrix_t -> filter[2][12] + *(j_1 + 13) * pMatrix_t -> filter[2][13] +
                                         *(j_1 + 14) * pMatrix_t -> filter[2][14] + *(j_1 + 15) * pMatrix_t -> filter[2][15] +
                                         *(j_1 + 16) * pMatrix_t -> filter[2][16] + *(j_1 + 17) * pMatrix_t -> filter[2][17];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;

        case 10:
            for(int i_num = 0;i_num < chn;i_num++,++i_1,++i_2)
            {
                std::vector<double>::iterator j_1 = i_1 -> second.begin();
                std::vector<double>::iterator j_2 = i_2 -> second.begin();

                if(i_1 -> first % 2 == 0 && i_1 != (ref_EX_map.end()--))//行标为偶数并且不是最后一个map的一次卷积
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_1 + 12) * pMatrix_t -> filter[2][12] + *(j_1 + 13) * pMatrix_t -> filter[2][13] +
                                         *(j_1 + 14) * pMatrix_t -> filter[2][14] + *(j_1 + 15) * pMatrix_t -> filter[2][15] +
                                         *(j_1 + 16) * pMatrix_t -> filter[2][16] + *(j_1 + 17) * pMatrix_t -> filter[2][17] +
                                         *(j_1 + 18) * pMatrix_t -> filter[2][18] + *(j_1 + 19) * pMatrix_t -> filter[2][19] +
                                         *(j_2 + 0) * pMatrix_t -> filter[3][0] + *(j_2 + 1) * pMatrix_t -> filter[3][1] +
                                         *(j_2 + 2) * pMatrix_t -> filter[3][2] + *(j_2 + 3) * pMatrix_t -> filter[3][3] +
                                         *(j_2 + 4) * pMatrix_t -> filter[3][4] + *(j_2 + 5) * pMatrix_t -> filter[3][5] +
                                         *(j_2 + 6) * pMatrix_t -> filter[3][6] + *(j_2 + 7) * pMatrix_t -> filter[3][7] +
                                         *(j_2 + 8) * pMatrix_t -> filter[3][8] + *(j_2 + 9) * pMatrix_t -> filter[3][9] +
                                         *(j_2 + 10) * pMatrix_t -> filter[3][10] + *(j_2 + 11) * pMatrix_t -> filter[3][11] +
                                         *(j_2 + 12) * pMatrix_t -> filter[3][12] + *(j_2 + 13) * pMatrix_t -> filter[3][13] +
                                         *(j_2 + 14) * pMatrix_t -> filter[3][14] + *(j_2 + 15) * pMatrix_t -> filter[3][15] +
                                         *(j_2 + 16) * pMatrix_t -> filter[3][16] + *(j_2 + 17) * pMatrix_t -> filter[3][17] +
                                         *(j_2 + 18) * pMatrix_t -> filter[3][18] + *(j_2 + 19) * pMatrix_t -> filter[3][19];         
                    }
                    ++i_1,++i_2;               
                }else if(i_1 -> first % 2 != 0 && i_1 != (ref_EX_map.end()--)){//行标为奇数并且不是最后一个map的一次卷积
                    
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1,++j_2)
                    {
                       ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[3][0] + *(j_1 + 1) * pMatrix_t -> filter[3][1] +
                                        *(j_1 + 2) * pMatrix_t -> filter[3][2] + *(j_1 + 3) * pMatrix_t -> filter[3][3] +
                                        *(j_1 + 4) * pMatrix_t -> filter[3][4] + *(j_1 + 5) * pMatrix_t -> filter[3][5] +
                                        *(j_1 + 6) * pMatrix_t -> filter[3][6] + *(j_1 + 7) * pMatrix_t -> filter[3][7] +
                                        *(j_1 + 8) * pMatrix_t -> filter[3][8] + *(j_1 + 9) * pMatrix_t -> filter[3][9] +
                                        *(j_1 + 10) * pMatrix_t -> filter[3][10] + *(j_1 + 11) * pMatrix_t -> filter[3][11] +
                                        *(j_1 + 12) * pMatrix_t -> filter[3][12] + *(j_1 + 13) * pMatrix_t -> filter[3][13] +
                                        *(j_1 + 14) * pMatrix_t -> filter[3][14] + *(j_1 + 15) * pMatrix_t -> filter[3][15] +
                                        *(j_1 + 16) * pMatrix_t -> filter[3][16] + *(j_1 + 17) * pMatrix_t -> filter[3][17] +
                                        *(j_1 + 18) * pMatrix_t -> filter[3][18] + *(j_1 + 19) * pMatrix_t -> filter[3][19];
                    }                    
                }else if(i_1 -> first % 2 == 0 && i_1 == (ref_EX_map.end()--))//行标为偶数并且是最后一个map的卷积 
                {
                    for(int j_int = 0; j_int < ref_CONLen;++j_int,++j_1)
                    {
                        ref_CON[j_int] = *(j_1 + 0) * pMatrix_t -> filter[2][0] + *(j_1 + 1) * pMatrix_t -> filter[2][1] +
                                         *(j_1 + 2) * pMatrix_t -> filter[2][2] + *(j_1 + 3) * pMatrix_t -> filter[2][3] +
                                         *(j_1 + 4) * pMatrix_t -> filter[2][4] + *(j_1 + 5) * pMatrix_t -> filter[2][5] +
                                         *(j_1 + 6) * pMatrix_t -> filter[2][6] + *(j_1 + 7) * pMatrix_t -> filter[2][7] +
                                         *(j_1 + 8) * pMatrix_t -> filter[2][8] + *(j_1 + 9) * pMatrix_t -> filter[2][9] +
                                         *(j_1 + 10) * pMatrix_t -> filter[2][10] + *(j_1 + 11) * pMatrix_t -> filter[2][11] +
                                         *(j_1 + 12) * pMatrix_t -> filter[2][12] + *(j_1 + 13) * pMatrix_t -> filter[2][13] +
                                         *(j_1 + 14) * pMatrix_t -> filter[2][14] + *(j_1 + 15) * pMatrix_t -> filter[2][15] +
                                         *(j_1 + 16) * pMatrix_t -> filter[2][16] + *(j_1 + 17) * pMatrix_t -> filter[2][17] +
                                         *(j_1 + 18) * pMatrix_t -> filter[2][18] + *(j_1 + 19) * pMatrix_t -> filter[2][19];
                    }       
                    ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON)); 
                    break;             
                }               
                ref_CON_map.insert(make_pair(i_1 -> first / 2,ref_CON));                 
            }
        break;     
    }
    //选择信号

    std::map<int,std::vector<double> >::iterator imp = ref_CON_map.begin();    
    ref_signal_map.erase(ref_signal_map.begin(),ref_signal_map.end());
    for(;imp!= ref_CON_map.end();++imp)
    {   
        ref_signal.erase(ref_signal.begin(),ref_signal.end());
        std::vector<double>::iterator ircm = imp -> second.begin();
        for(int i = 0; i < filterLen -1; ++i)
        {
            ircm++;
        }
        for(int i = 0; i < ref_dataLen; ++i,++ircm)
        {
            ref_signal.push_back(*ircm);
        }
        ref_signal_map.insert(make_pair(imp->first,ref_signal));
    }
    choose_row = ref_signal_map.size();
    REF_Len = ref_signal_map.begin() -> second .size();//重构信号长度重新定义

    if(DBN_N - 1 == n)
    {
        std::map<int,std::vector<double> >::iterator irsm = ref_signal_map.begin();
        std::vector<double>::iterator ivd = irsm -> second.begin();
        std::ofstream WD(DIR); 
        for(; ivd != irsm -> second.end(); ++ivd)
        {          
            WD << std::setprecision(16) << *ivd << " ";          
        }
    }
}
//延拓函数
void RTF::Reverse_transform_EX()
{    
    ref_EX = std::vector<double>(ref_EXLen);   
    //ref_EX_map.erase(ref_EX_map.begin(),ref_EX_map.end());       
    std::map<int,std::vector<double> >::iterator imus = ref_USam_map.begin(); //ref_signal_map的迭代器   

    for(std::vector<double>::size_type i = 0; i < static_cast<std::vector<double>::size_type>(choose_row); ++i,++imus)                             //全部延拓
    {       
        std::vector<double>::iterator idus = imus -> second.begin();
        for(std::vector<double>::size_type j = 0; j < static_cast<std::vector<double>::size_type>(filterLen - 1); ++j) //对称延拓
        {
            ref_EX[j]                               = *(idus + filterLen - 2 - j);//ref_USam[filterLen - 2 - j];
            ref_EX[filterLen + ref_USamLen - 1 + j] = *(idus + ref_USamLen - 1 - j);//ref_USam[ref_dataLen - 1 - j];
        }
        for(int j = 0; j != ref_USamLen; ++j)
        {
            ref_EX[filterLen - 1 + j] = *(idus + j);//ref_USam[j];
        }
        ref_EX_map.insert(make_pair(imus -> first, ref_EX));
    }
    return ;
}
//重构上采样函数
void RTF::Reverse_transform_USam()
{
    ref_USam = std::vector<double>(ref_USamLen);        
    ref_signal = std::vector<double>(REF_Len); 
    //ref_USam_map.erase(ref_USam_map.begin(),ref_USam_map.end());
    std::map<int,std::vector<double> >::iterator im = ref_signal_map.begin(); //ref_signal_map的迭代器 

    for(std::vector<double>::size_type i = 0; i < static_cast<std::vector<double>::size_type>(choose_row); ++i,++im)     //重构上采样赋值Alpha
    {
        ref_signal = im -> second;
        for(std::vector<double>::size_type j = 0; j < static_cast<std::vector<double>::size_type>(REF_Len); ++j)
        {
            ref_USam[2 * j] = 0; //ref_USam 行向量代码
            ref_USam[2 * j + 1] = ref_signal[j];
        }  
        ref_USam[ref_USamLen - 1] = 0;
        ref_USam_map.insert(make_pair(im -> first , ref_USam));  //制作usam map           
    }
    return;
}

void RTF::Reverse_transform_print(int n)                            //打印重构结果
{//0打印原始数据,1打印上采样数据,2是打印延拓数据,3是打印卷积数据,4是打印重构数据
    switch(n)
    {
        case 0:
            std::cout << "打印原始数据;"<< std::endl;
            std::cout << "是打印原始数据的宽度;"<< ref_signal_map.size() << std::endl;
            std::cout << "是打印原始数据的长度;"<< ref_signal_map.begin() -> second.size() << std::endl;
            std::cout << "是打印原始行标;"<< ref_signal_map.begin() -> first << std::endl;
            for(std::map<int,std::vector<double> >::iterator irsm = ref_signal_map.begin(); irsm != ref_signal_map.end();++irsm)
            {
                for(std::vector<double>::iterator irs = irsm -> second.begin();irs != irsm -> second.end(); ++irs)
                {   
                    std::cout << *irs << " ";
                }
                std::cout << std::endl;
            }
        break;

        case 1:
            std::cout << "1打印上采样数据;"<< std::endl;
            std::cout << "3是打印上采样数据的宽度;"<< ref_USam_map.size() << std::endl;
            std::cout << "3是打印上采样数据的长度;"<< ref_USam_map.begin() -> second.size() << std::endl;
            std::cout << "3是打印上采样行标;"<< ref_USam_map.begin() -> first << std::endl;
        /*    for(std::map<int,std::vector<double> >::iterator irsm = ref_USam_map.begin(); irsm != ref_USam_map.end();++irsm)
            {
                for(std::vector<double>::iterator irs = irsm -> second.begin();irs != irsm -> second.end(); ++irs)
                {   
                    std::cout << *irs << " ";
                }
                std::cout << std::endl;
            }
        */
       for(std::map<int,std::vector<double> >::iterator imus = ref_USam_map.begin();imus != ref_USam_map.end();++imus)
       {
           std::cout << imus -> first << std::endl;
       }
        break;

        case 2:
            std::cout << "2是打印延拓数据;"<< std::endl;
            std::cout << "3是打印延拓数据的宽度;"<< ref_EX_map.size() << std::endl;
            std::cout << "3是打印延拓数据的长度;"<< ref_EX_map.begin() -> second.size() << std::endl;
            std::cout << "3是打印延拓行标;"<< ref_EX_map.begin() -> first << std::endl;
        /*
            for(std::map<int,std::vector<double> >::iterator irsm = ref_EX_map.begin(); irsm != ref_EX_map.end();++irsm)
            {
                for(std::vector<double>::iterator irs = irsm -> second.begin();irs != irsm -> second.end(); ++irs)
                {   
                    std::cout << *irs << " ";
                }
                std::cout << std::endl;
            }
        */    
            for(std::map<int,std::vector<double> >::iterator imus = ref_EX_map.begin();imus != ref_EX_map.end();++imus)
            {
                std::cout << imus -> first << std::endl;
            }
        break;
        
        case 3:
            std::cout << "3是打印卷积数据;"<< std::endl;
            std::cout << "3是打印卷积数据的宽度;"<< ref_CON_map.size() << std::endl;
            std::cout << "3是打印卷积数据的长度;"<< ref_CON_map.begin() -> second.size() << std::endl;
            std::cout << "3是打印卷积行标;"<< ref_CON_map.begin() -> first << std::endl;
            for(std::map<int,std::vector<double> >::iterator irsm = ref_CON_map.begin(); irsm != ref_CON_map.end();++irsm)
            {
                for(std::vector<double>::iterator irs = irsm -> second.begin();irs != irsm -> second.end(); ++irs)
                {   
                    std::cout << *irs << " ";
                }
                std::cout << std::endl;
            }
        break;

        case 4:
            std::cout << "4是打印重构数据;"<< std::endl;
            for(std::map<int,std::vector<double> >::iterator irsm = ref_signal_map.begin(); irsm != ref_signal_map.end();++irsm)
            {
                for(std::vector<double>::iterator irs = irsm -> second.begin();irs != irsm -> second.end(); ++irs)
                {   
                    std::cout << *irs << " ";
                }
                std::cout << std::endl;
            }
        break;
    }

}

void RTF::Reverse_transform_Matrix(pMAT pMatrix)            //接受系数矩阵
{
    pMatrix_t = pMatrix; 
}
void RTF::Reverse_transform_REF_DATA(pREF pREF_data)         //接收分解的数据
{
    pREF_data_t = pREF_data;
}