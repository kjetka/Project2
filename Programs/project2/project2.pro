TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    non_interact.cpp

HEADERS += \
    non_interact.h

#INCLUDEPATH += D:\armadillo-8.100.1/include

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas

QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS -= -O2
