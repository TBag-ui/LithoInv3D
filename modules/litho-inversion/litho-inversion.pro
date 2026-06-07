QT -= gui
CONFIG += c++17 staticlib
TARGET = litho-inversion
TEMPLATE = lib

INCLUDEPATH += include \
    ../litho-core/include \
    ../litho-surface/include \
    ../litho-model/include \
    ../litho-forward/include \
    ../litho-em/include \
    ../litho-io/include \
    ../litho-regularization/include \
    ../../vendor/eigen

msvc: QMAKE_CXXFLAGS += /std:c++17
gcc: QMAKE_CXXFLAGS += -std=c++17

HEADERS += \
    include/litho_invert/inversion/constraint_handler.h \
    include/litho_invert/inversion/objective.h \
    include/litho_invert/inversion/joint_objective.h \
    include/litho_invert/inversion/optimizer.h \
    include/litho_invert/inversion/lbfgsb_optimizer.h \
    include/litho_invert/inversion/gncg_optimizer.h \
    include/litho_invert/inversion/runner.h

SOURCES += \
    src/inversion/constraint_handler.cpp \
    src/inversion/objective.cpp \
    src/inversion/joint_objective.cpp \
    src/inversion/lbfgsb_optimizer.cpp \
    src/inversion/gncg_optimizer.cpp \
    src/inversion/runner.cpp \
    src/inversion/ini_config_wiring.cpp
