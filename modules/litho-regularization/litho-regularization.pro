QT -= gui
CONFIG += c++17 staticlib
TARGET = litho-regularization
TEMPLATE = lib

INCLUDEPATH += include ../litho-core/include ../litho-surface/include ../litho-model/include ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

HEADERS += \
    include/litho_invert/regularization/regularization.h \
    include/litho_invert/regularization/smoothness.h \
    include/litho_invert/regularization/reference_model.h

SOURCES += \
    src/regularization/smoothness.cpp \
    src/regularization/reference_model.cpp
