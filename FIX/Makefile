TARGET = fixFD
CC = clang
CFLAGS = -O -I$(INCLUDE) 
INCLUDE = .

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).o $(LIBS)

clean:
	-rm $(TARGET).o

v:
	vi $(TARGET).c
