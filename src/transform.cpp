#include "head.h"

TF::TF(int dbn,int dbn_n,int row,int dataLen,int m,int n,int mode):DBN(dbn),DBN_N(dbn_n),data_row(row),data_line(dataLen),file_m(m),file_n(n),mode(mode)
{
    dec_dataLen = data_line;
    filterLen = 2 * DBN;
    dec_EXLen = dec_dataLen + 2 * filterLen - 2;
    dec_CONLen = dec_dataLen + filterLen - 1;
    dec_DSamLen = dec_CONLen / 2;
}
TF::~TF()
{


}
void TF::TransForm_mkdir_file()
{
    DIR = "./TF_DATE/s";
    std::string FILE_NAME = boost::lexical_cast<std::string>(file_m);
    DIR = DIR + FILE_NAME;
    DIR = "mkdir -p " + DIR;
    system(DIR.c_str());//创建Sn总文件夹
    std::string  num= boost::lexical_cast<std::string>(file_n);
    DIR = DIR + "/s"+ FILE_NAME+ "_" + num;
    system(DIR.c_str());//创建Sn_n总文件夹
}  
void TF::TransForm_mkdir_txt(int n)
{
    DIR = "./TF_DATE/s";
    std::string FILE_NAME = boost::lexical_cast<std::string>(file_m);
    std::string  num= boost::lexical_cast<std::string>(file_n);
    std::string  NEMBER_NAME= boost::lexical_cast<std::string>(n);
    DIR = DIR + FILE_NAME;// /home/lcx/code/REF_DATE/S1
    DIR = DIR + "/s"+ FILE_NAME+ "_" + num + "/s" + FILE_NAME + "_" + num + "_" + NEMBER_NAME + ".txt"; //  /s1/s1_1/s1_1_n.txt    
}
//dbn为分解层数,这个函数为小波分解主函数，用来调用其他子函数
void TF::transform(pRWD pRWdate,pMAT pMatrix)
{    
    TransForm_receive(pRWdate);     //获取原始数据指针
    TransForm_Matrix(pMatrix);      //获取矩阵系数
    TransForm_mkdir_file();         //创建文件夹
    
    for(int i = 0; i < data_row; ++i)
    {
        signal_row = std::vector<std::vector<double> >(1, std::vector<double>(data_line)); 
        TransForm_row(i); 
        TransForm_mkdir_txt(i + 1);                    //读取某一方数据 n = 0,代表获取第一行数据
        for(int j = 0; j < DBN_N; ++j)//分解DBN_N层
        {
            TransForm_EX(j);            //分解延拓
            TransForm_CON(j);           //分解卷积
            TransForm_DSam(j);          //分解下采样 
            //更新系数  
            dec_dataLen = DEC_Len;
            dec_EXLen = dec_dataLen + 2 * filterLen - 2;
            dec_CONLen = dec_dataLen + filterLen - 1;
            dec_DSamLen = dec_CONLen / 2;  
            
            if(filterLen >= dec_DSamLen)
            {
                std::cout << "滤波器长度大于信号长度！"<<std::endl;
                exit(-1);
            }            
               
        }
        //重新赋值
        dec_dataLen = data_line;
        filterLen = 2 * DBN;
        dec_EXLen = dec_dataLen + 2 * filterLen - 2;
        dec_CONLen = dec_dataLen + filterLen - 1;
        dec_DSamLen = dec_CONLen / 2;
    }    
    return ;
} 

/*
上采样函数得出分解结果
*/
void TF::TransForm_DSam(int n)
{
    signal_row = std::vector<std::vector<double> >(power(2,n + 1), std::vector<double>(dec_DSamLen));  
  
    if(n == DBN_N -1)
    {
        pDEC_t = new DEC;
        pDEC_t -> DEC_date = std::vector<std::vector<double> >(power(2,n + 1), std::vector<double>(dec_DSamLen));
       
        for(int i = 0; i < power(2,n); ++i)        //下采样
        {
            for(std::vector<double>::size_type j = 1; j <= static_cast<std::vector<double>::size_type>(dec_DSamLen); ++j)
            {
                signal_row[2 * i][j - 1]             = dec_CON[2 * i][2 * j - 1];
                signal_row[2 * i + 1][j - 1]         = dec_CON[2 * i +1][2 * j - 1];
                pDEC_t -> DEC_date[2 * i][j - 1]     = dec_CON[2 * i][2 * j - 1];
                pDEC_t -> DEC_date[2 * i + 1][j -1 ] = dec_CON[2 * i + 1][2 * j - 1];
               
            }
        }
        std::ofstream WD(DIR); 
        for(uint i = 0; i < signal_row.size(); ++i)
        {
            for(uint j = 0; j < signal_row[0].size(); ++j)
            {
                WD << std::setprecision(16) << signal_row[i][j] << " ";
            }
            WD << std::endl;
        }
               
        DEC_Len = dec_DSamLen;
    }else
    {
        for(int i = 0; i < power(2,n); ++i)        //下采样
        {
            for(std::vector<double>::size_type j = 1; j <= static_cast<std::vector<double>::size_type>(dec_DSamLen); ++j)
            {
                signal_row[2 * i][j - 1] = dec_CON[2 * i][2 * j - 1];
                signal_row[2 * i + 1][j - 1] = dec_CON[2 * i + 1][2 * j - 1];                
            }
        }
        DEC_Len = dec_DSamLen;
    }  

    return;
}

/*
    n代表延拓的层数，和分解的层数同步
*/
void TF::TransForm_CON(int n)//卷积函数
{
    dec_CON = std::vector<std::vector<double> >(power(2,n + 1), std::vector<double>(dec_CONLen));
    switch(DBN)
    {
        case 1:
            for(int i = 0; i < power(2,n); ++i)//卷积
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =  dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1];
                    dec_CON[2 * i + 1][j] =  dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1];
                }
            }
        break;
        
        case 2:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     = dec_EX[i][j]     * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] + 
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3]; 
                    dec_CON[2 * i + 1][j] = dec_EX[i][j]     * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] + 
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3]; 

                }                
            }
        break;

        case 3:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0;j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =     dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[0][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[0][5];
                    dec_CON[2 * i + 1][j] =     dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[1][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[1][5];            

                }

            }
        break;

        case 4:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =     dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[0][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[0][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[0][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[0][7] ;
                    dec_CON[2 * i + 1][j] =     dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[1][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[1][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[1][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[1][7] ; 
                }                           
            }
        break;

        case 5:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =      dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[0][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[0][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[0][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[0][7] +
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[0][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[0][9] ;
                    dec_CON[2 * i + 1][j] =     dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[1][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[1][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[1][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[1][7] + 
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[1][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[1][9] ;
                }                           
            }
        break;

        case 6:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =      dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[0][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[0][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[0][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[0][7] +
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[0][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[0][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[0][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[0][11] ;
                    dec_CON[2 * i + 1][j] =     dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[1][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[1][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[1][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[1][7] + 
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[1][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[1][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[1][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[1][11] ;
                }                           
            }
        break;

        case 7:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =      dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[0][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[0][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[0][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[0][7] +
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[0][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[0][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[0][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[0][11] +
                                            dec_EX[i][j + 12] * pMatrix_t -> filter[0][12] + dec_EX[i][j + 13] * pMatrix_t -> filter[0][13] ;
                    dec_CON[2 * i + 1][j] =     dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[1][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[1][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[1][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[1][7] + 
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[1][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[1][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[1][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[1][11] +
                                            dec_EX[i][j + 12] * pMatrix_t -> filter[1][12] + dec_EX[i][j + 13] * pMatrix_t -> filter[1][13] ;
                }                           
            }
        break;

        case 8:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =      dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[0][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[0][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[0][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[0][7] +
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[0][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[0][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[0][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[0][11] +
                                            dec_EX[i][j + 12] * pMatrix_t -> filter[0][12] + dec_EX[i][j + 13] * pMatrix_t -> filter[0][13] +
                                            dec_EX[i][j + 14] * pMatrix_t -> filter[0][14] + dec_EX[i][j + 15] * pMatrix_t -> filter[0][15] ;
                    dec_CON[2 * i + 1][j] =     dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[1][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[1][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[1][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[1][7] + 
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[1][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[1][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[1][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[1][11] +
                                            dec_EX[i][j + 12] * pMatrix_t -> filter[1][12] + dec_EX[i][j + 13] * pMatrix_t -> filter[1][13] +
                                            dec_EX[i][j + 14] * pMatrix_t -> filter[1][14] + dec_EX[i][j + 15] * pMatrix_t -> filter[1][15] ;
                }                           
            }
        break;

        case 9:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =      dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[0][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[0][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[0][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[0][7] +
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[0][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[0][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[0][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[0][11] +
                                            dec_EX[i][j + 12] * pMatrix_t -> filter[0][12] + dec_EX[i][j + 13] * pMatrix_t -> filter[0][13] +
                                            dec_EX[i][j + 14] * pMatrix_t -> filter[0][14] + dec_EX[i][j + 15] * pMatrix_t -> filter[0][15] +
                                            dec_EX[i][j + 16] * pMatrix_t -> filter[0][16] + dec_EX[i][j + 17] * pMatrix_t -> filter[0][17] ;
                    dec_CON[2 * i + 1][j] =     dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[1][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[1][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[1][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[1][7] + 
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[1][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[1][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[1][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[1][11] +
                                            dec_EX[i][j + 12] * pMatrix_t -> filter[1][12] + dec_EX[i][j + 13] * pMatrix_t -> filter[1][13] +
                                            dec_EX[i][j + 14] * pMatrix_t -> filter[1][14] + dec_EX[i][j + 15] * pMatrix_t -> filter[1][15] +
                                            dec_EX[i][j + 16] * pMatrix_t -> filter[1][16] + dec_EX[i][j + 17] * pMatrix_t -> filter[1][17] ;
                }                           
            }
        break;

        case 10:
            for(int i = 0; i < power(2,n); ++i)
            {
                for(int j = 0; j < dec_CONLen; ++j)
                {
                    dec_CON[2 * i][j]     =      dec_EX[i][j] * pMatrix_t -> filter[0][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[0][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[0][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[0][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[0][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[0][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[0][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[0][7] +
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[0][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[0][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[0][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[0][11] +
                                            dec_EX[i][j + 12] * pMatrix_t -> filter[0][12] + dec_EX[i][j + 13] * pMatrix_t -> filter[0][13] +
                                            dec_EX[i][j + 14] * pMatrix_t -> filter[0][14] + dec_EX[i][j + 15] * pMatrix_t -> filter[0][15] +
                                            dec_EX[i][j + 16] * pMatrix_t -> filter[0][16] + dec_EX[i][j + 17] * pMatrix_t -> filter[0][17] +
                                            dec_EX[i][j + 18] * pMatrix_t -> filter[0][18] + dec_EX[i][j + 19] * pMatrix_t -> filter[0][19] ;
                    dec_CON[2 * i + 1][j] =     dec_EX[i][j] * pMatrix_t -> filter[1][0] + dec_EX[i][j + 1] * pMatrix_t -> filter[1][1] +
                                            dec_EX[i][j + 2] * pMatrix_t -> filter[1][2] + dec_EX[i][j + 3] * pMatrix_t -> filter[1][3] + 
                                            dec_EX[i][j + 4] * pMatrix_t -> filter[1][4] + dec_EX[i][j + 5] * pMatrix_t -> filter[1][5] +
                                            dec_EX[i][j + 6] * pMatrix_t -> filter[1][6] + dec_EX[i][j + 7] * pMatrix_t -> filter[1][7] + 
                                            dec_EX[i][j + 8] * pMatrix_t -> filter[1][8] + dec_EX[i][j + 9] * pMatrix_t -> filter[1][9] +
                                            dec_EX[i][j + 10] * pMatrix_t -> filter[1][10] + dec_EX[i][j + 11] * pMatrix_t -> filter[1][11] +
                                            dec_EX[i][j + 12] * pMatrix_t -> filter[1][12] + dec_EX[i][j + 13] * pMatrix_t -> filter[1][13] +
                                            dec_EX[i][j + 14] * pMatrix_t -> filter[1][14] + dec_EX[i][j + 15] * pMatrix_t -> filter[1][15] +
                                            dec_EX[i][j + 16] * pMatrix_t -> filter[1][16] + dec_EX[i][j + 17] * pMatrix_t -> filter[1][17] +
                                            dec_EX[i][j + 18] * pMatrix_t -> filter[1][18] + dec_EX[i][j + 19] * pMatrix_t -> filter[1][19] ;
                }                           
            }
        break;
    }   
    return;
}

/*
    n代表延拓的层数，和分解的层数同步
*/
void TF::TransForm_EX(int n)//分解延拓函数,n为分解的
{  
    dec_EX = std::vector<std::vector<double> >(power(2,n), std::vector<double>(dec_EXLen));

    for(std::vector<double>::size_type i = 0; i < static_cast<std::vector<double>::size_type>(power(2,n)); ++i)                             //全部延拓
    {
        for(std::vector<double>::size_type j = 0; j < static_cast<std::vector<double>::size_type>(filterLen - 1); ++j) //对称延拓
        {
            dec_EX[i][j]            = signal_row[i][filterLen - 2 - j];
            dec_EX[i][filterLen + dec_dataLen - 1 + j] = signal_row[i][dec_dataLen - 1 - j];
        }
    }
    for(std::vector<double>::size_type i = 0; i != signal_row.size(); ++i)
    {
        for(std::vector<double>::size_type j = 0; j != signal_row[i].size(); ++j)
        {
            dec_EX[i][filterLen - 1 + j] = signal_row[i][j];
        }
    }
    return;
}

void TF::TransForm_row(int n)//获取单一通道数据,n读取n行的数据
{      
    for(int i = 0; i < data_line; ++i)
    {        
        signal_row[0][i] = dec_pRWdate -> Raw_Data[n][i];
    }

    return;
}

/*
n = 0,打印矩阵原始数据，n = 1，打印读取的行向量。n = 2,打印延拓函数n = 3打印卷积函数,n = 4打印采样函数
*/
void TF::TransForm_print(int n)
{
    switch (n)
    {
        case 0:
            std::cout << "Print_TransForm:"<<  dec_pRWdate -> Raw_Data.size() << std::endl;
            for(std::vector<double>::size_type i = 0; i != dec_pRWdate -> Raw_Data.size(); ++i)
            {
                for(std::vector<double>::size_type j = 0; j != dec_pRWdate -> Raw_Data[i].size(); ++j)
                {
                    //std::cout << "dec_pRWdate -> Raw_Data[" << i <<"]["<< j << "] = " << dec_pRWdate -> Raw_Data[i][j] << "   ";
                    std::cout << dec_pRWdate -> Raw_Data[i][j] << " ";
                }
                std::cout<< std::endl;
            }
        break;

        case 1:            
            for(std::vector<double>::size_type i = 0; i != signal_row.size(); ++i)
            {
                for(std::vector<double>::size_type j = 0; j != signal_row[i].size(); ++j)
                {
                    //std::cout << "signal_row[" << i << "]["<< j << "] = " << signal_row[i][j] << "    "; 
                    std::cout << signal_row[i][j] << " ";
                }                
            }
            std::cout << std::endl;
        break;

        case 2:            
            for(std::vector<double>::size_type i = 0; i != dec_EX.size(); ++i)
            {
                for(std::vector<double>::size_type j = 0; j != dec_EX[i].size(); ++j)
                {
                    //std::cout << "dec_EX[" << i << "]["<< j << "] = " << dec_EX[i][j] << "    "; 
                    std::cout << dec_EX[i][j] << " ";
                }                
            }
            std::cout << std::endl;
        break;

        case 3:            
            for(std::vector<double>::size_type i = 0; i != dec_CON.size(); ++i)
            {
                for(std::vector<double>::size_type j = 0; j != dec_CON[i].size(); ++j)
                {
                    //std::cout << "dec_CON[" << i << "]["<< j << "] = " << dec_CON[i][j] << "    "; 
                    std::cout << dec_CON[i][j] << " ";
                }                
            }
            std::cout << std::endl;
        break;

        case 4:   
        
            for(std::vector<double>::size_type i = 0; i != signal_row.size(); ++i)
            {
                for(std::vector<double>::size_type j = 0; j !=signal_row[i].size(); ++j)
                {
                    //std::cout << "dec_DSam[" << i << "]["<< j << "] = " << signal_row[i][j] << "    "; 
                    std::cout << signal_row[i][j] << " ";
                }      
                std::cout << std::endl;          
            }            
        break;        
    }

    return ;
}
void TF::TransForm_receive(pRWD pRWdate)//获取读到的数据
{
    dec_pRWdate = pRWdate;
}
void TF::TransForm_Matrix(pMAT pMatrix)
{
    pMatrix_t = pMatrix;
}
