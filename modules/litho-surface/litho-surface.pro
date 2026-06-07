QT -= gui
CONFIG += c++17 staticlib
TARGET = litho-surface
TEMPLATE = lib

INCLUDEPATH += include ../litho-core/include ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

HEADERS += \
    include/litho_invert/surface/surface_mesh.h \
    include/litho_invert/surface/boundary_loop.h

SOURCES += \
    src/surface/surface_mesh.cpp \
    src/surface/boundary_loop.cpp
