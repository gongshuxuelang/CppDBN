.PYONY:clean
cc = g++
CFLAGS = -Wall -g -I ./inc -I. -std=c++11  #头文件
LIBO = -lpthread -lm	#链接
BIN = main
SUBDIR = $(shell ls -d */)  
ROOTSRC = $(wildcard *.cpp)     
ROOTOBJ = $(ROOTSRC:%.c=%.o)
SUBSRC = $(shell find $(SUBDIR) -name '*.cpp')
SUBOBJ = $(SUBSRC:%.c=%.o)
 
$(BIN):$(ROOTOBJ) $(SUBOBJ)
		$(cc) $(CFLAGS) -o $(BIN) $(ROOTOBJ)  $(SUBOBJ) $(LIBO)
 $(SUBOBJ).c.o:
	 $(cc) $(CFLAGS) -c $< -o $@
 
clean:
	  rm -f $(BIN) $(ROOTOBJ) $(SUBOBJ) 



