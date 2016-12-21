bPath:=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/

exec:=main

ALL: main.o
	g++ -g  -llapack main.o -o $(exec)


main.o: main.cpp
	g++ -g -I $(bPath) main.cpp -c
run: 
	./$(exec)
clean:
	rm -rf *.o $(exec)
