QT -= gui
CONFIG += c++17 staticlib
TARGET = litho-em
TEMPLATE = lib

INCLUDEPATH += include ../litho-core/include ../litho-surface/include ../litho-model/include ../litho-forward/include ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

HEADERS += \
    include/litho_invert/em/em_data.h \
    include/litho_invert/em/em_solver.h \
    include/litho_invert/em/em_forward_model.h \
    include/litho_invert/em/em_active_forward.h \
    include/litho_invert/em/em_mt_forward.h

SOURCES += \
    src/em/em_forward_model.cpp \
    src/em/em_solver.cpp \
    src/em/em_active_forward.cpp \
    src/em/em_mt_forward.cpp
