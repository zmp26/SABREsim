#ifndef CONSOLECOLORIZER_H
#define CONSOLECOLORIZER_H

#include <string>

class ConsoleColorizer{
public:
	static void PrintRed(const std::string& text);
	static void PrintGreen(const std::string& text);
	static void PrintBlue(const std::string& text);
	static void PrintYellow(const std::string& text);
	static void PrintOrange(const std::string& text);
	static void PrintWhite(const std::string& text);

	static void PrintCode(const std::string& text, const std::string& code);
	static void PrintRGB(const std::string& text, int r, int g, int b);
};

#endif