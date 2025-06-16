// config.cpp
#include "config.hpp"
#include <filesystem>
#include <iostream>
#include <cstdlib>

namespace fs = std::filesystem;

const std::vector<std::string> colorPalette = {
    "lightblue", "lightcoral", "lightgreen", "khaki",
    "plum", "orange", "turquoise", "wheat",
    "thistle", "gold", "pink", "lightsalmon",
    "lightgray", "darkseagreen", "lavender"
};

Config& Config::getInstance() {
    static Config instance;
    return instance;
}


void Config::setOutputDir(const std::string& dir) {
    // if (fs::exists(dir)) {
    //     std::cerr << "Error: Output directory \"" << dir << "\" already exists.\n";
    //     std::exit(EXIT_FAILURE);
    // }

    try {
        fs::create_directories(dir);
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error: Failed to create output directory: " << e.what() << "\n";
        std::exit(EXIT_FAILURE);
    }

    outputDir = dir;

    // Open log file
    std::string logPath = outputDir + "/log.txt";
    logFile.open(logPath, std::ios::out | std::ios::app);
    if (!logFile.is_open()) {
        std::cerr << "Error: Could not open log file: " << logPath << "\n";
        std::exit(EXIT_FAILURE);
    }
}


std::string Config::getOutputDir() const {
    return outputDir;
}


void Config::log(const std::string& message) {
    if (logFile.is_open()) {
        logFile << message;
    } else {
        std::cerr << "Log file not open." << std::endl;
    }
}
