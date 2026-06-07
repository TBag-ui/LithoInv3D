QT -= gui
CONFIG += c++17 console
TARGET = litho-model-tests
TEMPLATE = app

INCLUDEPATH += \
    ../include \
    ../../litho-core/include \
    ../../litho-surface/include \
    ../../../vendor/eigen \
    ../../../vendor

msvc: QMAKE_CXXFLAGS += /std:c++17 /EHsc
gcc: QMAKE_CXXFLAGS += -std=c++17

LIBS += \
    -L../../../build/release \
    -llitho-core \
    -llitho-surface \
    -llitho-model

SOURCES += \
    test_no_gap.cpp \
    test_global_dofs.cpp \
    test_fix_exterior.cpp \
    test_apply_bounds.cpp \
    ../../../vendor/catch2/catch_amalgamated.cpp
