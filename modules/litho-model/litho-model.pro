QT -= gui
CONFIG += c++17 staticlib
TARGET = litho-model
TEMPLATE = lib

INCLUDEPATH += include ../litho-core/include ../litho-surface/include ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

HEADERS += \
    include/litho_invert/litho/litho_group.h \
    include/litho_invert/litho/lithology_model.h \
    include/litho_invert/inversion/inversion_result.h

SOURCES += \
    src/litho/litho_group.cpp \
    src/litho/lithology_model.cpp
