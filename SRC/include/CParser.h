
#ifndef __CPARSER_H_
#define __CPARSER_H_

    
#include <typeinfo>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <fstream>

using namespace std;


//void exitWithError(const std::string &error) {
//    std::cout << error;
//    std::cin.ignore();
//    std::cin.get();
//    //exit(EXIT_FAILURE);
//}

// -----------------------------------------------------------------------------
// Conversion of std::string to primitive types (int/float/double/...), 
// and vice-versa
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Conversion of std::string to primitive types (int/float/double/...), 
// and vice-versa
// -----------------------------------------------------------------------------
class Convert {

    public:

        // Convert T, which should be a primitive, to a std::string
        template <typename T>
        static std::string T_to_string(T const &val) {
            std::ostringstream ostr;
            ostr << val;
            return ostr.str();
        }

        // Convert a std::string to T.
        // I have specialized string_to_T function, for std::string. 
        // Why? Well, take a look at this: if (!(istr >> returnVal))
        // if function parameter val would be a string containing 
        // whitespace, like: toyota corolla -> then string_to_T will 
        // return only "toyota", since istringstream will stop \
        // extracting  at the first whitespace.
        template <typename T>
        static T string_to_T(std::string const &val) {
            std::istringstream istr(val);
            T returnVal;
            if (!(istr >> returnVal))
                cout << "CFG : Not a valid type: " + (std::string)typeid(T).name() + " received!\n" << endl;
    //            exitWithError("CFG: Not a valid " + (std::string)typeid(T).name() + " received!\n");
            return returnVal;

        }
      //  template <>
      //  static std::string string_to_T(std::string const &val) {
      //      return val;
      //  }

}; // end Convert --------------------------------------------------------------






// -----------------------------------------------------------------------------
// Contains functions needed to parse the configuration file. 
// -----------------------------------------------------------------------------
class ConfigFile {

    private:

        // store pairs key-value
        std::map<std::string, std::string> contents;
        // config name
        std::string fName;

        // remove comments from config file
        // function resturns false if a non-space character is found
        // if so it removes it removes everything from the semicolon
        // to the end of the line.
        void removeComment(std::string &line) const {
            if (line.find(';') != line.npos) line.erase(line.find(';'));
        }
        // removes free lines (sometimes, after remove a comment, the 
        // line only contains whitespaces...)
        // Basically, the function returns false if a non-space character 
        // was found, true otherwise. The function is "const" because it 
        // does not alter any class member variables.
        bool onlyWhitespace(const std::string &line) const {
            return (line.find_first_not_of(' ') == line.npos);
        }
        // checks if an individual line has the correct structure of 
        // a config file (key = value)
        // the function accepts as parameter a std::string, which is an 
        // individual line (with removed comment), from the config file.
        bool validLine(const std::string &line) const {
            std::string temp = line;
            // The .erase() simply removes every character starting from 
            // position 0 -> first non-tab or non-space character. After 
            // removal, if the first character is '=', it means that we 
            // do not have a key. Something like this: "  =someValue"
            temp.erase(0, temp.find_first_not_of("\t "));
            if (temp[0] == '=') return false;
            // The for loop loops starting from the position of the '=', 
            // until the end of the line. If a non-space character was found, 
            // then we have a key value. If the "if" never executes, the 
            // function returns false, because the key doesn't have a value. An example in which the function also returns false
            for (size_t i = temp.find('=') + 1; i < temp.length(); i++)
                if (temp[i] != ' ') return true;
            return false;
        }

        // extracts the key from the pair of key = value.
        // sepPos represents the position of the '='
        // Example: (   car = ford ) -> The value of key will be: "   car" (there are 
        // three whitespaces in front of car, but DIC code tags won't allow 
        // whitespace). Why? Because that .substr() creates a substring 
        // starting with the character at position 0, and finishes with the 
        // character from the position of '=' - 1. Then, everything from 
        // the first space or tab character, is removed.
        void extractKey(std::string &key, size_t const &sepPos, const std::string &line) const    {
            key = line.substr(0, sepPos);
            if (key.find('\t') != line.npos || key.find(' ') != line.npos)
                key.erase(key.find_first_of("\t "));
        }

        // Extracts the key
        // Example (car = toyota corolla) value will be assigned 
        // "toyota corolla". First of all, .substr() creates a substring 
        // starting from positon of '=' + 1, to the end of the line. 
        // Then, value.erase(0, value.find_first_not_of("\t ")); removes the 
        // leading whitespace, and value.erase(value.find_last_not_of("\t ") + 1); 
        // removes everything starting with the position of the last non-tab or non-space character.
        void extractValue(std::string &value, size_t const &sepPos, const std::string &line) const {
            value = line.substr(sepPos + 1);
            value.erase(0, value.find_first_not_of("\t "));
            value.erase(value.find_last_not_of("\t ") + 1);
        }


        /// Following functions calls the above funcions ///


        // To store in map member "contents" a pair key-value
        void extractContents(const std::string &line)  {
            std::string temp = line;
            temp.erase(0, temp.find_first_not_of("\t "));
            size_t sepPos = temp.find('=');
            std::string key, value;
            extractKey(key, sepPos, temp);
            extractValue(value, sepPos, temp);
            if (!keyExists(key))
                contents.insert(std::pair<std::string, std::string>(key, value));
            else {
            
                //exitWithError("CFG: Can only have unique key names!\n");
                cout << "CFG: Can only have unique key names!\n" << endl;
                std::cout << "Sometimes after using \"Align in vim\" does not work" << std::endl;
            }
        }

        // Calls the above function
        void parseLine(const std::string &line, size_t const lineNo) {
            if (line.find('=') == line.npos)
                //exitWithError("CFG: Couldn't find separator on line: " + Convert::T_to_string(lineNo) + "\n");
                cout << "CFG: Couldn't find separator on line: " + Convert::T_to_string(lineNo) + "\n" << endl;
            if (!validLine(line))
                //exitWithError("CFG: Bad format for line: " + Convert::T_to_string(lineNo) + "\n");
                cout << "CFG: Couldn't find separator on line: " + Convert::T_to_string(lineNo) + "\n" << endl;
            extractContents(line);
        }


    private:
        // opens the configuration file, and extracts & parses it's contents
        void ExtractKeys()  {
            std::ifstream file;
            file.open(fName.c_str());
            if (!file)
                //exitWithError("CFG: File " + fName + " couldn't be found!\n");
                cout << "CFG: File " + fName + " couldn't be found!\n" << endl;
            std::string line;
            size_t lineNo = 0;
            // loop keeps extracting lines, until EOF is found
            while (std::getline(file, line)) {
                lineNo++;
                std::string temp = line;
                if (temp.empty())
                    continue;
                removeComment(temp);
                if (onlyWhitespace(temp))
                    continue;
                parseLine(temp, lineNo);
            }
            file.close();
        }


    public:

        ConfigFile(const std::string &fName);
        
        // function which keys if a specific key exists in the configuration 
        // file. Since the pair of key-value is extracted in our map, all we 
        // have to do is to use std::map::find function to look for the key:
        bool keyExists(const std::string &key) const {
            return contents.find(key) != contents.end();
        }

        // retrieves the value of a specific key
        // returns a default value (operator()() of ValueType), if the key 
        // couldn't be found. Otherwise, it will return the converted value 
        // from string to ValueType, of the key. 
        template <typename ValueType>
        ValueType getValueOfKey(const std::string &key, ValueType const &defaultValue = ValueType()) const {
            if (!keyExists(key))
                return defaultValue;
            return Convert::string_to_T<ValueType>(contents.find(key)->second);
        }

};

#endif
