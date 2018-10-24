QT       -= core gui
CONFIG   += c++11
QMAKE_CXXFLAGS += -std=c++11

TARGET = brainnets
TEMPLATE = lib
VERSION = 1.0.0
DEFINES += BRAINNETS_LIBRARY

QMAKE_LN_SHLIB = :
CONFIG += unversioned_libname unversioned_soname skip_target_version_ext
win32:TARGET_EXT = .dll

SOURCES += brainnets.cpp \
    ts.cpp \
    networks.cpp \
    networks_mst.cpp \
    distances.cpp

HEADERS += brainnets.h\
        brainnets_global.h \
    ts.h \
    networks.h \
    networks_mst.h \
    sample_data.h \
    distances.h \
    ../extra/eigen_extensions.h \

DISTFILES += \
    README.md

INSTALLS += target
INCLUDEPATH += $$PWD/../../third-party/
INCLUDEPATH += $$PWD/../../third-party/boost/
INCLUDEPATH += $$PWD/../../third-party/eigen3/

unix {
    LIBS += -lgomp
}

CONFIG(debug, debug|release) {
    DESTDIR = ../../builds/debug/
    OBJECTS_DIR = ../../builds/debug/brainnets-tmp/
    MOC_DIR = ../../builds/debug/brainnets-tmp/
    RCC_DIR = ../../builds/debug/brainnets-tmp/
    UI_DIR = ../../builds/debug/brainnets-tmp/
} else {
    DESTDIR = ../../builds/release/
    OBJECTS_DIR = ../../builds/release/brainnets-tmp/
    MOC_DIR = ../../builds/release/brainnets-tmp/
    RCC_DIR = ../../builds/release/brainnets-tmp/
    UI_DIR = ../../builds/release/brainnets-tmp/
}
