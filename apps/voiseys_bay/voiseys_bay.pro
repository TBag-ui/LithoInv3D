QT -= gui
CONFIG += c++17 console
TARGET = csv_invert_fixed
TEMPLATE = app

INCLUDEPATH += \
    ../../modules/litho-core/include \
    ../../modules/litho-surface/include \
    ../../modules/litho-model/include \
    ../../modules/litho-forward/include \
    ../../modules/litho-em/include \
    ../../modules/litho-io/include \
    ../../modules/litho-regularization/include \
    ../../modules/litho-inversion/include \
    ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

LIBS += \
    -L../../build/release \
    -llitho-core \
    -llitho-surface \
    -llitho-model \
    -llitho-forward \
    -llitho-em \
    -llitho-io \
    -llitho-regularization \
    -llitho-inversion

SOURCES += main.cpp generate_synthetic.cpp
HEADERS += generate_synthetic.h
