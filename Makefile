# TARGET := main
TARGET := main

run: build prepare
	./$(TARGET).o;

build:
	g++ -o $(TARGET).o $(TARGET).cpp -O2 -lm -lpthread -I/usr/X11R6/include -L/usr/X11R6/lib -lm -lpthread -lX11;

prepare:
	mkdir -p result;

clean:
	@echo "Cleaning..."
	rm -rf $(TARGET).o
	rm -rf result