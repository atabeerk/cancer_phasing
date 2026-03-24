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

    void setMaxHopMultiplier(int value);
    int getMaxHopMultiplier() const;

private:
    Config() = default;
    std::string outputDir;
    std::ofstream logFile;
    int maxHopMultiplier = 4;
};

extern const std::vector<std::string> colorPalette;


#endif // CONFIG_H
