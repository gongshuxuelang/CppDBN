src存储的是源文件，inc存储的是头文件。makefile是脚本文件。

本次备份是完成小波一次分解和合成实现的通用代码。
构造函数有两个参数需要赋值，第一个参数是进行dbn，例如进行db2，第一个参数输入2，进行db4，第二个参数输入4.目前支持 n = 1~10 。第二个参数是数据的位数。比如信号的个数是8，第二个参数就输入8，信号个数是8192，第二个参数的个数就是8192.注意信号必须是2的n次幂的形式。
Copyright (C) 2019 Li Chunxue

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
