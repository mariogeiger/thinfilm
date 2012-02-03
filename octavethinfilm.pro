QT       += core
QT       -= gui

TARGET = multilayer-thinfilm-octave
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

SOURCES += octave.cpp
HEADERS += thinfilm.hpp

target.path = /usr/local/bin
INSTALLS += target
