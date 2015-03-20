CXX = g++
CXXFLAGS = -ansi -pedantic -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -std=c++0x -DLINUX -fPIE -m64
TARGET = qttest
BUILD_PATH = build
SRC_PATH= src
INCLUDE_PATH = -I include -isystem /opt/Qt/5.4/gcc_64/include
LIBS_PATH =  -L /opt/Qt/5.4/gcc_64/lib
LIBS= -lQt5Widgets -lQt5Test -lQt5Gui -lQt5Core -Wl,-rpath=/opt/Qt/Tools/QtCreator/lib/qtcreator

UIC = uic
UI_FORMS = ui_mainwindow.h

OBJS = testwindow.o  main.o
F_OBJS = $(addprefix $(BUILD_PATH)/, $(OBJS))
MOCS = moc_testwindow.cpp
MOCS_OBJS = moc_testwindow.o
F_MOCS_OBJS = $(addprefix $(BUILD_PATH)/, $(MOCS_OBJS))


vpath %.ui junk/untitled/
vpath %.cpp src
vpath %.h include
vpath %.cpp build

.PONY: all clean run

all: $(UI_FORMS) $(MOCS) $(MOCS_OBJS) $(OBJS)
	$(CXX)  $(F_OBJS) $(F_MOCS_OBJS) -o $(BUILD_PATH)/$(TARGET) $(LIBS_PATH) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS)  $(INCLUDE_PATH) -c $< -o $(BUILD_PATH)/$@

moc_%.cpp: %.h
	moc $< -o  $(BUILD_PATH)/$@

clean:
	rm -f ui/*.ui
	rm -f include/ui_*.h
	rm -f build/*.o
	rm -f build/$(TARGET)

run: clean all
	./$(BUILD_PATH)/$(TARGET) 

ui_%.h: %.ui
	cp $< ui/$(notdir $<)
	$(UIC) $< -o include/$@