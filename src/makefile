CC = gcc
CFLAGS = -O2 -Wall --std=c11 -I./include -lm
COMOBJ = ./lib/config.c 
TARGET1 = even
TARGET2 = prop

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(COMOBJ) ./lib/$(TARGET1).c
	$(CC) -o $@ $^ $(CFLAGS)

prop: $(COMOBJ) ./lib/$(TARGET2).c
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	$(RM) $(TARGET1) $(TARGET2)
