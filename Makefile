
PROF =
WARN = -Wall -Wextra -Wuninitialized -Wlogical-op 
DFLAGS = -g -DS_DEBUG -ftrapv
OFLAGS = 
CFLAGS = $(WARN) $(DFLAGS) $(OFLAGS) $(PROF) -pipe 
IFLAGS =
LFLAGS = 
LIBS = -lstdc++ -lm

DEP_SRC = *.cc

SOURCE = *.cc *.h

OBJECTS = main.o configfile.o msdatafile.o msumsoptions.o runanalysis.o statheader.o stats_single.o

TARGET = msums
DIRECTORY = msums

all : $(TARGET)

include Makefile.common

include Makefile.dep

install : $(TARGET)
	cp $(TARGET) $(HOME)/bin

install_strip : install
	strip $(HOME)/bin/$(TARGET)

release :
	$(MAKE) DFLAGS="" OFLAGS="-O3 $(ARCH)" && strip $(TARGET)

clean :
	rm -f *.o *.moc

all_clean : clean
	rm -f $(TARGET)
