QT -= gui
CONFIG += c++17 staticlib
TARGET = litho-forward
TEMPLATE = lib

INCLUDEPATH += include ../litho-core/include ../litho-surface/include ../litho-model/include ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

HEADERS += \
    include/litho_invert/forward/forward_model.h \
    include/litho_invert/forward/gravity_forward.h \
    include/litho_invert/forward/magnetic_forward.h

SOURCES += \
    src/forward/gravity_forward.cpp \
    src/forward/magnetic_forward.cpp
