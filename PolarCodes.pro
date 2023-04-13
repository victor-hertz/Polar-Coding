QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -Ofast

SOURCES += \
        AWGNChannel.cpp \
        AbstractChannel.cpp \
        AbstractPolarCode.cpp \
        BSChannel.cpp \
        PolarCode.cpp \
        Utils.cpp \
        main.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    AWGNChannel.h \
    AbstractChannel.h \
    AbstractPolarCode.h \
    BSChannel.h \
    PolarCode.h \
    Utils.h
