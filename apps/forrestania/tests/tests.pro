QT -= gui
CONFIG += c++17 console
TARGET = forrestania-tests
TEMPLATE = app

INCLUDEPATH += \
    .. \
    ../../../modules/litho-core/include \
    ../../../modules/litho-surface/include \
    ../../../modules/litho-model/include \
    ../../../modules/litho-forward/include \
    ../../../modules/litho-em/include \
    ../../../modules/litho-io/include \
    ../../../modules/litho-regularization/include \
    ../../../modules/litho-inversion/include \
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
    test_forrestania_integration.cpp \
    ../generate_synthetic.cpp \
    ../cluster_loader.cpp \
    ../../../vendor/catch2/catch_amalgamated.cpp
