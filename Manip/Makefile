BASE = ../../rai-fork
BASE2 = ../../rai-manip
NAME   = $(shell basename `pwd`)
OUTPUT = lib$(NAME).so

DEPEND = Core Algo Kin Gui Geo KOMO Optim Control

EIGEN = 1
OMPL = 1

#CPATHS += /usr/include/ompl
CPATHS += /usr/include/eigen3

#LPATHS += /usr/local/lib/libompl.so
#LIBS += /usr/local/lib/libompl.so

#DEPEND_UBUNTU += libompl-dev
#LIBS += -lompl -lboost_regex

#LIBS += -fopenmp -unroll-loops
#LIBS += -unroll-loops

LIBS += -lspdlog -lfmt  -g 
LPATHS += -L/usr/lib/x86_64-linux-gnu/libspdlog.so

SRCS = $(shell find . -maxdepth 2 -name '*.cpp' )
OBJS = $(SRCS:%.cpp=%.o)

include $(BASE)/build/generic.mk
