#   make        -- compiles your project into program.exe
#   make clean  -- removes compiled item
#
# All .cpp flles are included.

CC = mpic++
SRCS = $(wildcard *.cpp)
HDRS = $(wildcard *.h)
OBJS = $(SRCS:.cpp=.o)
DIRS = $(subst /, ,$(CURDIR))
PROJ = mpi_optics

APP = $(PROJ)
CFLAGS= -c -O2 -w 
LDFLAGS=
LIBS= -lstdc++ -O2

all: $(APP)

$(APP): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(APP) $(LIBS)

%.o: %.cpp $(HDRS) $(MF)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(APP)

