#ifndef BRAINNETS_GLOBAL_H
#define BRAINNETS_GLOBAL_H

#if defined(BRAINNETS_LIBRARY)
    #if defined(WIN32)
        #define BRAINNETSSHARED_EXPORT __declspec(dllexport)
    #else
        #define BRAINNETSSHARED_EXPORT __attribute__ ((visibility ("default")))
    #endif
#else
    #if defined(WIN32)
        #define BRAINNETSSHARED_EXPORT __declspec(dllimport)
    #else
        #define BRAINNETSSHARED_EXPORT __attribute__ ((visibility ("default")))
    #endif
#endif

namespace brainnets {

#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define VERSION_PATCH 0
#define VERSION_EXTRA "alpha"

}

#endif // BRAINNETS_GLOBAL_H
