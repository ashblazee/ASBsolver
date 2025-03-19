#include <iostream>
#include <iomanip>  // For setting width and percentages
#include <chrono>

// Function to show progress bar
void showProgressBar(long long current, long long total, long double value) {
    int barWidth = 50;  // Width of the loading bar
    float progress = (float)current / total;
    int pos = barWidth * progress;
    std::cout << value << ' ';
    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";     // Completed portion
        else if (i == pos) std::cout << ">"; // Current point
        else std::cout << " ";             // Remaining portion
    }
    std::cout << "] " << std::setw(3) << int(progress * 100.0) << "%\r";
    std::cout.flush();
}

