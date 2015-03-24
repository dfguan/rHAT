CC=g++
CFLAGS=-O3 -c -Wall
DFLAGS=-g -Wall
LDFLAGS= -lz -lm
ECPPSOURCES=main.cpp aligner.cpp graph.cpp  readfl.cpp hash.cpp form.cpp 
CSOURCES=ksw.c 
EOBJECTS=$(ECPPSOURCES:.cpp=.o) $(CSOURCES:.c=.o)
EXECUTABLE=rHAT-mapper
HCPPSOURCES=whash.cpp genhash.cpp formh.cpp
HOBJECTS=$(HCPPSOURCES:.cpp=.o) $(CSOURCES:.c=.o)
HEXECUTABLE=rHAT-indexer

all: $(HCPPSOURCES) $(ECPPSOURCES) $(CSOURCES) $(EXECUTABLE) $(HEXECUTABLE)
	
$(EXECUTABLE): $(EOBJECTS)
	$(CC)  $(EOBJECTS) -o $@ $(LDFLAGS)

$(HEXECUTABLE): $(HOBJECTS)
	$(CC)  $(HOBJECTS) readfl.o -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
.c.o:
	$(CC) $(CFLAGS) $< -o $@
debugE:
	$(CC) $(DFLAGS) $(ECPPSOURCES) $(CSOURCES) $(LDFLAGS)
debugH:
	$(CC) $(DFLAGS) $(HCPPSOURCES) $(CSOURCES) readfl.cpp $(LDFLAGS)
printT:
	$(CC) $(DFLAGS) -pg $(ECPPSOURCES) $(CSOURCES) $(LDFLAGS)
clean:
	rm -rf *.o $(EXECUTABLE) $(HEXECUTABLE) a.out

 
