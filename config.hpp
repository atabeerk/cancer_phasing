// config.h
#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <fstream>
#include <vector>

class Config {
public:
    static Config& getInstance();

    void setOutputDir(const std::string& dir);
    std::string getOutputDir() const;
    void log(const std::string& message);

private:
    Config() = default;
    std::string outputDir;
    std::ofstream logFile;
};

extern const std::vector<std::string> colorPalette;


#endif // CONFIG_H
