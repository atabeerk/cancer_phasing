#include "pipeline.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    try {
        PipelineOptions options;
        const auto result = parsePipelineArguments(argc, argv, options);
        if (result == ArgumentResult::Help) return 0;
        if (result == ArgumentResult::Error) return 1;
        return runChromosomePipeline(options);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}
