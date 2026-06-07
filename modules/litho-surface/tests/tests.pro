QT -= gui
CONFIG += c++17 console
TARGET = litho-surface-tests
TEMPLATE = app

INCLUDEPATH += \
    ../include \
    ../../litho-core/include \
    ../../../vendor/eigen \
    ../../../vendor

msvc: QMAKE_CXXFLAGS += /std:c++17 /EHsc
gcc: QMAKE_CXXFLAGS += -std=c++17

LIBS += \
    -L../../../build/release \
    -llitho-core \
    -llitho-surface

SOURCES += \
    test_vertex_freedom.cpp \
    test_dof_system.cpp \
    test_face_normals.cpp \
    ../../../vendor/catch2/catch_amalgamated.cpp
