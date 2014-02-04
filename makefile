_OBJS =	bilayer.o tools.o sample.o g2model.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = g++
ODIR = obj
SDIR = src
IDIR = inc

CFLAGS = -Wall -O3 -ffast-math -g
CPATH = -I/usr/include/tcl8.4 -I./$(IDIR)
CLINK = -ltcl8.4 -ltk8.4 -lBLT24

all: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(CPATH) $(CLINK) -o sdp
	
$(ODIR)/%.o: $(SDIR)/%.cpp $(IDIR)/%.h
	$(CC) -c $(CPATH) $(CFLAGS) -o $@ $< 
	
clean:
	rm -f $(ODIR)/*.o sdp
