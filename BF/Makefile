TARGET = bfFD
CC = clang
CFLAGS = -O -I$(INCLUDE) #### -DDEBUG
INCLUDE = .

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).o $(LIBS)

clean:
	-rm $(TARGET).o $(TARGET)

v:
	vi $(TARGET).c
