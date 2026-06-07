QT -= gui
CONFIG += c++17 staticlib
TARGET = litho-io
TEMPLATE = lib

INCLUDEPATH += include ../litho-core/include ../litho-surface/include ../litho-model/include ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

HEADERS += \
    include/litho_invert/io/io.h \
    include/litho_invert/io/gravity_loader.h \
    include/litho_invert/io/surface_loader.h \
    include/litho_invert/io/constraint_loader.h \
    include/litho_invert/io/litho_config_loader.h \
    include/litho_invert/io/ini_config.h \
    include/litho_invert/io/dem_loader.h \
    include/litho_invert/io/exporters.h \
    include/litho_invert/io/magnetic_loader.h \
    include/litho_invert/plot/svg_plot.h

SOURCES += \
    src/io/gravity_loader.cpp \
    src/io/surface_loader.cpp \
    src/io/constraint_loader.cpp \
    src/io/litho_config_loader.cpp \
    src/io/ini_config.cpp \
    src/io/dem_loader.cpp \
    src/io/exporters.cpp \
    src/io/magnetic_loader.cpp
