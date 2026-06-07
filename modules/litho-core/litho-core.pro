QT -= gui
CONFIG += c++17 staticlib
TARGET = litho-core
TEMPLATE = lib

INCLUDEPATH += include ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

HEADERS += \
    include/litho_invert/litho_invert.h \
    include/litho_invert/core/common.h \
    include/litho_invert/core/geometry.h \
    include/litho_invert/core/stats.h

SOURCES += \
    src/core/common.cpp \
    src/core/geometry.cpp \
    src/core/stats.cpp
