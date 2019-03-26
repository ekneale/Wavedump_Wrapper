# Check for command goals other than clean.
ifeq ($(MAKECMDGOALS),clean)
	CLEAN_RESULT:=$(shell rm -f $(NS)/src/*.o $(NS)/src/*.d $(NS)/src/*~ *~ *.o *.d)
	ERROR_RESULT:=$(error "Finished cleaning project.")
endif

ifeq ($(MAKECMDGOALS),)
	ERROR_RESULT:=$(error "No executable name provided: make (name)")
endif

EXTRAOBJS = $(MAKECMDGOALS).o

# g++ is the compiler.
CC = g++

# Directories for all include files.
IFOLDER = -I$(NS)/src `root-config --cflags`

# Directories for all libraries.
LFOLDER = 

# Check operating system.
UNAME = $(shell uname)

# Libraries to link to.
#LIBS = -lgsl -lgslcblas `root-config --glibs`
LIBS =  `root-config --glibs`
# Option flags for all build modes.
CCOPT = -g -O2 -Wall -pedantic -Wno-long-long

# Name of executable.
TARGET = $(MAKECMDGOALS)

# List of object files.
OBJS = $(patsubst %.cpp,%.o,$(wildcard $(NS)/src/*.cpp)) $(EXTRAOBJS)

# List of dependencies.
DEPS = $(patsubst %.o,%.d,$(OBJS))

$(MAKECMDGOALS): $(OBJS) $(DEPS)
	$(CC) $(CCOPT) -o $(TARGET) $(LFOLDER) $(OBJS) $(LIBS)

%.o: %.cpp
	$(CC) $(CCOPT) $(IFOLDER) -o $@ -c $<

%.d: %.cpp
	$(CC) $(CCOPT) $(IFOLDER) -MM -MT '$(subst Q,$(patsubst %.cpp,%,$<),Q.o Q.d)' $< > $@

clean:
	rm -f $(NS)/src/*.o $(NS)/src/*.d $(NS)/src/*~ *~

# If the command is not make clean, generate dependencies.
ifneq (clean,$(findstring $(MAKECMDGOALS),clean))
-include $(DEPS)
endif

