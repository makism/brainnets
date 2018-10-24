#ifndef BRAINNETS_H
#define BRAINNETS_H

// STD/STL
#include <string>

// BrainNets
#include "brainnets_global.h"
#include "networks.h"
#include "networks_mst.h"
#include "ts.h"
#include "sample_data.h"

namespace brainnets {

std::string get_full_version() {
    return std::to_string(VERSION_MAJOR) + "." +
            std::to_string(VERSION_MINOR) + "." +
            std::to_string(VERSION_PATCH) + "-" +
            VERSION_EXTRA;

}

}

#endif // BRAINNETS_H
