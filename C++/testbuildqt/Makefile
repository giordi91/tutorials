#compiler varaiables 
CXX = g++
CXXFLAGS = -ansi -pedantic -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -std=c++0x -DLINUX -fPIE -m64
#needed libs
LIBS_PATH =  -L /opt/Qt/5.4/gcc_64/lib
LIBS= -lQt5Widgets -lQt5Test -lQt5Gui -lQt5Core -Wl,-rpath=/opt/Qt/Tools/QtCreator/lib/qtcreator
#includes variables
INCLUDE_PATH = -I include -isystem /opt/Qt/5.4/gcc_64/include

#build variables, name of the application and build path
TARGET = qttest
BUILD_PATH = build
SRC_PATH= src

#UI compiler variables
UIC = uic
UI_FORMS = ui_mainwindow.h
#Mocs varaibles
MOCS = moc_testwindow.cpp
MOCS_OBJS = moc_testwindow.o
F_MOCS_OBJS = $(addprefix $(BUILD_PATH)/, $(MOCS_OBJS))
#objs to be generated and linked
OBJS = testwindow.o  main.o
F_OBJS = $(addprefix $(BUILD_PATH)/, $(OBJS))

#setup of the search paths
vpath %.ui junk/untitled/
vpath %.cpp src
vpath %.h include
vpath %.cpp build

#pony targets
.PONY: all clean run

#in chager of running needed generic compiling rules and linking
all: $(UI_FORMS) $(MOCS) $(MOCS_OBJS) $(OBJS)
	$(CXX)  $(F_OBJS) $(F_MOCS_OBJS) -o $(BUILD_PATH)/$(TARGET) $(LIBS_PATH) $(LIBS)

#generic rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS)  $(INCLUDE_PATH) -c $< -o $(BUILD_PATH)/$@

moc_%.cpp: %.h
	moc $< -o  $(BUILD_PATH)/$@

ui_%.h: %.ui
	cp $< ui/$(notdir $<)
	$(UIC) $< -o include/$@

#clean the project
clean:
	rm -f ui/*.ui
	rm -f include/ui_*.h
	rm -f build/*.o
	rm -f build/$(TARGET)

#target for cleaning compiling and running the project
run: clean all
	./$(BUILD_PATH)/$(TARGET) 
