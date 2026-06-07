QT -= gui
CONFIG += c++17 console
TARGET = litho-io-tests
TEMPLATE = app

INCLUDEPATH += \
    ../include \
    ../../litho-core/include \
    ../../litho-surface/include \
    ../../litho-model/include \
    ../../../vendor/eigen \
    ../../../vendor

msvc: QMAKE_CXXFLAGS += /std:c++17 /EHsc
gcc: QMAKE_CXXFLAGS += -std=c++17

LIBS += \
    -L../../../build/release \
    -llitho-core \
    -llitho-surface \
    -llitho-model \
    -llitho-io

SOURCES += \
    test_ts_content.cpp \
    test_no_duplicate_writers.cpp \
    ../../../vendor/catch2/catch_amalgamated.cpp
