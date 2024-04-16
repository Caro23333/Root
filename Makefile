vpath %.cpp src
vpath %.o build
vpath %.h include

obj = testlib.o CP_RootAnother.o CP_RootSturmNewton.o VAS.o CP_Root.o testMain.o
# header = GameWindow.h GameSocket.h Resource.h ResourceManager.h

all: testMain.exe

.PHONY: all clean

testMain.exe: $(obj)
	g++ $(addprefix build/,$(obj)) -o bin/testMain.exe -L ./lib -lgmpxx -lgmp -Wall

testlib.o: testlib.cpp testlib.h CP_Root.h
	g++ -c src/testlib.cpp -o build/testlib.o -I ./include -Wall

CP_RootAnother.o: CP_RootAnother.cpp CP_Root.h
	g++ -c src/CP_RootAnother.cpp -o build/CP_RootAnother.o -I ./include -Wall

CP_RootSturmNewton.o: CP_RootSturmNewton.cpp CP_Root.h testlib.h
	g++ -c src/CP_RootSturmNewton.cpp -o build/CP_RootSturmNewton.o -I ./include -Wall

VAS.o: VAS.cpp VAS.h
	g++ -c src/VAS.cpp -o build/VAS.o -I ./include -Wall

CP_Root.o: CP_Root.cpp CP_Root.h testlib.h
	g++ -c src/CP_Root.cpp -o build/CP_Root.o -I ./include -Wall

testMain.o: test/testMain.cpp testlib.h CP_Root.h
	g++ -c test/testMain.cpp -o build/testMain.o -I ./include -Wall

clean:
	-del /f /q build\* bin\*
