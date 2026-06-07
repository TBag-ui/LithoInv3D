QT -= gui
CONFIG += c++17 console
TARGET = litho-inversion-tests
TEMPLATE = app

INCLUDEPATH += \
    ../include \
    ../../litho-core/include \
    ../../litho-surface/include \
    ../../litho-model/include \
    ../../litho-forward/include \
    ../../litho-em/include \
    ../../litho-io/include \
    ../../litho-regularization/include \
    ../../../vendor/eigen \
    ../../../vendor

msvc: QMAKE_CXXFLAGS += /std:c++17 /EHsc
gcc: QMAKE_CXXFLAGS += -std=c++17

LIBS += \
    -L../../../build/release \
    -llitho-core \
    -llitho-surface \
    -llitho-model \
    -llitho-forward \
    -llitho-em \
    -llitho-io \
    -llitho-regularization \
    -llitho-inversion

SOURCES += \
    test_simple_runner.cpp \
    test_vertex_movement.cpp \
    test_model_bounds.cpp \
    test_iteration_ts_export.cpp \
    ../../../vendor/catch2/catch_amalgamated.cpp
