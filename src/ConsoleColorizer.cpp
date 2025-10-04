#include <iostream>
#include "ConsoleColorizer.h"

void ConsoleColorizer::PrintCode(const std::string& text, const std::string& code){
	std::cout << "\033[" << code << "m" << text << "\033[0m";
}

void ConsoleColorizer::PrintRGB(const std::string& text, int r, int g, int b){
	std::cout << "\033[38;2;" << r << ";" << g << ";" << b << "m" << text << "\033[0m";
}

void ConsoleColorizer::PrintRed(const std::string& text){
	//std::cout << "\033[31m" << text << "\033[0m";
	//PrintCode(text, "31");
	PrintRGB(text, 255, 0, 0);
}

void ConsoleColorizer::PrintGreen(const std::string& text){
	//std::cout << "\033[32m" << text << "\033[0m";
	//PrintCode(text, "32");
	PrintRGB(text, 0, 255, 0);
}

void ConsoleColorizer::PrintBlue(const std::string& text){
	//std::cout << "\033[34m" << text << "\033[0m";
	//PrintCode(text, "34");
	PrintRGB(text, 0, 0, 255);
}

void ConsoleColorizer::PrintYellow(const std::string& text){
	//std::cout << "\033[33m" << text << "\033[0m";
	//PrintCode(text, "33");
	PrintRGB(text, 255, 255, 0);
}

void ConsoleColorizer::PrintOrange(const std::string& text){
	//orange is not a default color w/ associated code, so let's use PrintRGB to do this!
	PrintRGB(text,255,165,0);
}

void ConsoleColorizer::PrintWhite(const std::string& text){
	PrintRGB(text,255,255,255);
}

void ConsoleColorizer::PrintPurple(const std::string& text){
	PrintRGB(text,70,29,124);
}

void ConsoleColorizer::PrintGold(const std::string& text){
	PrintRGB(text,253,208,35);
}