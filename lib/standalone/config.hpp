// ======================================================================
// Atomistica - Interatomic potential library and molecular dynamics code
// https://github.com/Atomistica/atomistica
//
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// and others. See the AUTHORS file in the top-level Atomistica directory.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ======================================================================

#pragma once

#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace atomistica {

// Simple ptrdict-style configuration file parser.
//
// Parses files of the form:
//   Simulation {
//     key = "value";
//     Section {
//       key = "value";
//     };
//   };
class Config {
public:
    void parse_file(const std::string& filename) {
        std::ifstream f(filename);
        if (!f)
            throw std::runtime_error("Cannot open config file: " + filename);
        std::string content((std::istreambuf_iterator<char>(f)),
                             std::istreambuf_iterator<char>());
        parse_block(content, 0);
    }

    double get_double(const std::string& key, double default_val = 0.0) const {
        auto it = values_.find(key);
        if (it == values_.end()) return default_val;
        return std::stod(it->second);
    }

    int get_int(const std::string& key, int default_val = 0) const {
        auto it = values_.find(key);
        if (it == values_.end()) return default_val;
        return std::stoi(it->second);
    }

    std::string get_string(const std::string& key,
                           const std::string& default_val = "") const {
        auto it = values_.find(key);
        if (it == values_.end()) return default_val;
        return it->second;
    }

    bool has_section(const std::string& name) const {
        // case-insensitive match
        std::string lower_name = to_lower(name);
        for (auto& kv : sections_) {
            if (to_lower(kv.first) == lower_name) return true;
        }
        return false;
    }

    const Config& section(const std::string& name) const {
        std::string lower_name = to_lower(name);
        for (auto& kv : sections_) {
            if (to_lower(kv.first) == lower_name) return kv.second;
        }
        throw std::runtime_error("Section not found: " + name);
    }

private:
    std::map<std::string, std::string> values_;
    std::map<std::string, Config> sections_;

    static std::string to_lower(const std::string& s) {
        std::string r = s;
        for (auto& c : r) c = static_cast<char>(std::tolower(c));
        return r;
    }

    static std::string strip(const std::string& s) {
        size_t a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return "";
        size_t b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    }

    // Remove comments (# to end-of-line)
    static std::string remove_comments(const std::string& src) {
        std::string out;
        out.reserve(src.size());
        bool in_str = false;
        for (size_t i = 0; i < src.size(); ++i) {
            char c = src[i];
            if (c == '"') in_str = !in_str;
            if (!in_str && c == '#') {
                // skip to end of line
                while (i < src.size() && src[i] != '\n') ++i;
                if (i < src.size()) out += '\n';
            } else {
                out += c;
            }
        }
        return out;
    }

    // Parse content inside a { ... } block.
    // Returns the position just after the closing '}'.
    // Used recursively for nested sections.
    size_t parse_block(const std::string& src, size_t pos) {
        std::string clean = remove_comments(src);
        return parse_block_inner(clean, pos);
    }

    size_t parse_block_inner(const std::string& src, size_t pos) {
        // Skip optional outer section name and opening brace for top-level call
        // Find the Simulation block if present, otherwise parse as flat
        size_t n = src.size();

        // Skip whitespace
        auto skip_ws = [&](size_t p) {
            while (p < n && std::isspace(src[p])) ++p;
            return p;
        };

        // Read an identifier (letters, digits, underscore, slash, dot)
        auto read_ident = [&](size_t p) -> std::pair<std::string, size_t> {
            std::string id;
            while (p < n && (std::isalnum(src[p]) || src[p] == '_' ||
                              src[p] == '/' || src[p] == '.'))
            {
                id += src[p++];
            }
            return {id, p};
        };

        pos = skip_ws(pos);

        // If we're at top level, look for "Simulation {" wrapper
        if (pos == 0) {
            // Find "Simulation"
            auto [name, after_name] = read_ident(pos);
            if (to_lower(name) == "simulation") {
                pos = skip_ws(after_name);
                if (pos < n && src[pos] == '{') {
                    ++pos; // consume '{'
                    pos = parse_inner(src, pos);
                    return pos;
                }
            }
            // No Simulation wrapper — parse as flat
            pos = 0;
        }

        return parse_inner(src, pos);
    }

    // Parse key=value pairs and Section { } blocks until '}' or end.
    // Returns position after the closing '}'.
    size_t parse_inner(const std::string& src, size_t pos) {
        size_t n = src.size();

        auto skip_ws = [&](size_t p) {
            while (p < n && std::isspace(src[p])) ++p;
            return p;
        };

        auto read_ident = [&](size_t p) -> std::pair<std::string, size_t> {
            std::string id;
            while (p < n && (std::isalnum(src[p]) || src[p] == '_' ||
                              src[p] == '/' || src[p] == '.'))
            {
                id += src[p++];
            }
            return {id, p};
        };

        while (pos < n) {
            pos = skip_ws(pos);
            if (pos >= n) break;

            if (src[pos] == '}') {
                ++pos; // consume '}'
                // skip optional ';'
                pos = skip_ws(pos);
                if (pos < n && src[pos] == ';') ++pos;
                break;
            }

            // Read identifier
            auto [ident, after_ident] = read_ident(pos);
            if (ident.empty()) {
                ++pos; // skip unknown character
                continue;
            }
            pos = skip_ws(after_ident);

            if (pos < n && src[pos] == '{') {
                // It's a section
                ++pos; // consume '{'
                Config sub;
                pos = sub.parse_inner(src, pos);
                sections_[ident] = std::move(sub);
            } else if (pos < n && src[pos] == '=') {
                // key = "value";
                ++pos; // consume '='
                pos = skip_ws(pos);
                std::string value;
                if (pos < n && src[pos] == '"') {
                    ++pos; // consume opening '"'
                    while (pos < n && src[pos] != '"') value += src[pos++];
                    if (pos < n) ++pos; // consume closing '"'
                } else {
                    // unquoted value — read until ';'
                    while (pos < n && src[pos] != ';' && src[pos] != '\n')
                        value += src[pos++];
                    value = strip(value);
                }
                values_[ident] = value;
                // consume optional ';'
                pos = skip_ws(pos);
                if (pos < n && src[pos] == ';') ++pos;
            } else {
                // Unexpected — skip to next ';' or '{'
                while (pos < n && src[pos] != ';' && src[pos] != '{' &&
                       src[pos] != '}') ++pos;
                if (pos < n && src[pos] == ';') ++pos;
            }
        }
        return pos;
    }
};

} // namespace atomistica
