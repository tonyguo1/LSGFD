CXX=mpicxx
CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

INCS= -I. 

OBJS =		LSGFD.o 

LIBS =

TARGET =	LSGFD

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
