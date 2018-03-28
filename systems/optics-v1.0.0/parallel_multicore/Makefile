#   make        -- compiles your project into program.exe
#   make clean  -- removes compiled item
#   make handin -- creates a project Zip file for hand in
#
# All .cpp flles are included.

CC = cc
SRCS = $(wildcard *.cpp)
HDRS = $(wildcard *.h)
OBJS = $(SRCS:.cpp=.o)
DIRS = $(subst /, ,$(CURDIR))
#PROJ = $(word $(words $(DIRS)), $(DIRS))
PROJ = omp_optics

APP = $(PROJ)
#CFLAGS = -c -g -Wall -I ../boost_1_47_0/
CFLAGS= -c -fopenmp -O2 -w #-I/usr/common/usg/pnetcdf/1.3.0/gnu/47/include
#-I ../../../kdtree2.31/src-c++/
LDFLAGS=
#LIBS= -lnetcdf_c++ -lnetcdf -fopenmp -O2 
#LIBS= -lnetcdf -fopenmp -O2 
LIBS= -fopenmp -O2 #-L /usr/common/usg/pnetcdf/1.3.0/gnu/47/lib -lpnetcdf

all: $(APP)

$(APP): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(APP) $(LIBS)

%.o: %.cpp $(HDRS) $(MF)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(APP)

run:
	./$(APP) -d 10 -p 100 -m 2 -e 0.8 -pr 1 

handin:
	rm -f $(PROJ).zip
	zip $(PROJ).zip -q Makefile *.cpp *.h
	unzip -l $(PROJ).zip

