CC = gcc
CFLAGS = -Wall -Wextra -lm

TARGET = audio_plugin
SRC = src/audio_plugin.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
