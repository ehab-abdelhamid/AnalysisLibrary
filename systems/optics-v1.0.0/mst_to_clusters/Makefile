#   make        -- compiles your project into program.exe
#   make clean  -- removes compiled item
#   make handin -- creates a project Zip file for hand in
#
# All .cpp flles are included.

CC = g++
SRCS = $(wildcard *.cpp)
HDRS = $(wildcard *.h)
OBJS = $(SRCS:.cpp=.o)
DIRS = $(subst /, ,$(CURDIR))
PROJ = mst_to_clusters

APP = $(PROJ)
CFLAGS= -c -fopenmp -O2 -w 
LDFLAGS=
LIBS= -fopenmp -O2 

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

