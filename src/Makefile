#------------------------------------------------------------------------------#
# This makefile was generated by 'cbp2make' tool rev.147                       #
#------------------------------------------------------------------------------#


WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = 
CFLAGS = -Wall -fexceptions
RESINC = 
LIBDIR = 
LIB = 
LDFLAGS = 

INC_DEBUG = $(INC)
CFLAGS_DEBUG = $(CFLAGS) -g
RESINC_DEBUG = $(RESINC)
RCFLAGS_DEBUG = $(RCFLAGS)
LIBDIR_DEBUG = $(LIBDIR)
LIB_DEBUG = $(LIB)
LDFLAGS_DEBUG = $(LDFLAGS)
OBJDIR_DEBUG = obj/Debug
DEP_DEBUG = 
OUT_DEBUG = bin/Debug/nnmerge

INC_RELEASE = $(INC) -I/usr/local/opt/boost1.58/include/
CFLAGS_RELEASE = $(CFLAGS) -O2 -std=c++11
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -s -lboost_timer -lboost_system
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/nnmerge

OBJ_DEBUG = $(OBJDIR_DEBUG)/smerge.o $(OBJDIR_DEBUG)/lysearch.o $(OBJDIR_DEBUG)/main.o $(OBJDIR_DEBUG)/scriptparser.o $(OBJDIR_DEBUG)/timer.o $(OBJDIR_DEBUG)/vmath.o $(OBJDIR_DEBUG)/vstring.o $(OBJDIR_DEBUG)/abstractmerge.o $(OBJDIR_DEBUG)/cleaner.o $(OBJDIR_DEBUG)/distopt.o $(OBJDIR_DEBUG)/evaluator.o $(OBJDIR_DEBUG)/abstractgraph.o $(OBJDIR_DEBUG)/iodelegator.o $(OBJDIR_DEBUG)/jmerge.o $(OBJDIR_DEBUG)/kgraph.o $(OBJDIR_DEBUG)/lynndescent.o

OBJ_RELEASE = $(OBJDIR_RELEASE)/smerge.o $(OBJDIR_RELEASE)/lysearch.o $(OBJDIR_RELEASE)/main.o $(OBJDIR_RELEASE)/scriptparser.o $(OBJDIR_RELEASE)/timer.o $(OBJDIR_RELEASE)/vmath.o $(OBJDIR_RELEASE)/vstring.o $(OBJDIR_RELEASE)/abstractmerge.o $(OBJDIR_RELEASE)/cleaner.o $(OBJDIR_RELEASE)/distopt.o $(OBJDIR_RELEASE)/evaluator.o $(OBJDIR_RELEASE)/abstractgraph.o $(OBJDIR_RELEASE)/iodelegator.o $(OBJDIR_RELEASE)/jmerge.o $(OBJDIR_RELEASE)/kgraph.o $(OBJDIR_RELEASE)/lynndescent.o

all: debug release

clean: clean_debug clean_release

before_debug: 
	test -d bin/Debug || mkdir -p bin/Debug
	test -d $(OBJDIR_DEBUG) || mkdir -p $(OBJDIR_DEBUG)

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(LD) $(LIBDIR_DEBUG) -o $(OUT_DEBUG) $(OBJ_DEBUG)  $(LDFLAGS_DEBUG) $(LIB_DEBUG)

$(OBJDIR_DEBUG)/smerge.o: smerge.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c smerge.cpp -o $(OBJDIR_DEBUG)/smerge.o

$(OBJDIR_DEBUG)/lysearch.o: lysearch.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c lysearch.cpp -o $(OBJDIR_DEBUG)/lysearch.o

$(OBJDIR_DEBUG)/main.o: main.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c main.cpp -o $(OBJDIR_DEBUG)/main.o

$(OBJDIR_DEBUG)/scriptparser.o: scriptparser.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c scriptparser.cpp -o $(OBJDIR_DEBUG)/scriptparser.o

$(OBJDIR_DEBUG)/timer.o: timer.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c timer.cpp -o $(OBJDIR_DEBUG)/timer.o

$(OBJDIR_DEBUG)/vmath.o: vmath.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c vmath.cpp -o $(OBJDIR_DEBUG)/vmath.o

$(OBJDIR_DEBUG)/vstring.o: vstring.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c vstring.cpp -o $(OBJDIR_DEBUG)/vstring.o

$(OBJDIR_DEBUG)/abstractmerge.o: abstractmerge.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c abstractmerge.cpp -o $(OBJDIR_DEBUG)/abstractmerge.o

$(OBJDIR_DEBUG)/cleaner.o: cleaner.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c cleaner.cpp -o $(OBJDIR_DEBUG)/cleaner.o

$(OBJDIR_DEBUG)/distopt.o: distopt.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c distopt.cpp -o $(OBJDIR_DEBUG)/distopt.o

$(OBJDIR_DEBUG)/evaluator.o: evaluator.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c evaluator.cpp -o $(OBJDIR_DEBUG)/evaluator.o

$(OBJDIR_DEBUG)/abstractgraph.o: abstractgraph.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c abstractgraph.cpp -o $(OBJDIR_DEBUG)/abstractgraph.o

$(OBJDIR_DEBUG)/iodelegator.o: iodelegator.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c iodelegator.cpp -o $(OBJDIR_DEBUG)/iodelegator.o

$(OBJDIR_DEBUG)/jmerge.o: jmerge.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c jmerge.cpp -o $(OBJDIR_DEBUG)/jmerge.o

$(OBJDIR_DEBUG)/kgraph.o: kgraph.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c kgraph.cpp -o $(OBJDIR_DEBUG)/kgraph.o

$(OBJDIR_DEBUG)/lynndescent.o: lynndescent.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c lynndescent.cpp -o $(OBJDIR_DEBUG)/lynndescent.o

clean_debug: 
	rm -f $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf bin/Debug
	rm -rf $(OBJDIR_DEBUG)

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/smerge.o: smerge.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c smerge.cpp -o $(OBJDIR_RELEASE)/smerge.o

$(OBJDIR_RELEASE)/lysearch.o: lysearch.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c lysearch.cpp -o $(OBJDIR_RELEASE)/lysearch.o

$(OBJDIR_RELEASE)/main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o $(OBJDIR_RELEASE)/main.o

$(OBJDIR_RELEASE)/scriptparser.o: scriptparser.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c scriptparser.cpp -o $(OBJDIR_RELEASE)/scriptparser.o

$(OBJDIR_RELEASE)/timer.o: timer.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c timer.cpp -o $(OBJDIR_RELEASE)/timer.o

$(OBJDIR_RELEASE)/vmath.o: vmath.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c vmath.cpp -o $(OBJDIR_RELEASE)/vmath.o

$(OBJDIR_RELEASE)/vstring.o: vstring.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c vstring.cpp -o $(OBJDIR_RELEASE)/vstring.o

$(OBJDIR_RELEASE)/abstractmerge.o: abstractmerge.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c abstractmerge.cpp -o $(OBJDIR_RELEASE)/abstractmerge.o

$(OBJDIR_RELEASE)/cleaner.o: cleaner.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c cleaner.cpp -o $(OBJDIR_RELEASE)/cleaner.o

$(OBJDIR_RELEASE)/distopt.o: distopt.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c distopt.cpp -o $(OBJDIR_RELEASE)/distopt.o

$(OBJDIR_RELEASE)/evaluator.o: evaluator.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c evaluator.cpp -o $(OBJDIR_RELEASE)/evaluator.o

$(OBJDIR_RELEASE)/abstractgraph.o: abstractgraph.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c abstractgraph.cpp -o $(OBJDIR_RELEASE)/abstractgraph.o

$(OBJDIR_RELEASE)/iodelegator.o: iodelegator.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c iodelegator.cpp -o $(OBJDIR_RELEASE)/iodelegator.o

$(OBJDIR_RELEASE)/jmerge.o: jmerge.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c jmerge.cpp -o $(OBJDIR_RELEASE)/jmerge.o

$(OBJDIR_RELEASE)/kgraph.o: kgraph.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c kgraph.cpp -o $(OBJDIR_RELEASE)/kgraph.o

$(OBJDIR_RELEASE)/lynndescent.o: lynndescent.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c lynndescent.cpp -o $(OBJDIR_RELEASE)/lynndescent.o

clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release

