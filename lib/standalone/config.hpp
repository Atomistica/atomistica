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
#include <utility>
#include <vector>

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
//
// Changes from the original version:
//   - sections_ is a std::vector<std::pair<std::string,Config>> so that
//     (a) declaration order is preserved for registry-driven dispatch and
//     (b) multiple sections with the same name (e.g., two OutputXYZ blocks)
//         are both retained.
//   - all_sections() exposes the ordered section list.
//   - get_or<T>() provides a uniform templated accessor.
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

    // -----------------------------------------------------------------------
    // Value accessors (key=value pairs in this section)
    // -----------------------------------------------------------------------

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

    bool get_bool(const std::string& key, bool default_val = false) const {
        auto it = values_.find(key);
        if (it == values_.end()) return default_val;
        std::string s = to_lower(it->second);
        return s == "true" || s == "yes" || s == "1";
    }

    // Uniform template accessor. Specialisations below.
    template<typename T>
    T get_or(const std::string& key, T default_val) const;

    bool has_key(const std::string& key) const {
        return values_.find(key) != values_.end();
    }

    // -----------------------------------------------------------------------
    // Sub-section accessors
    // -----------------------------------------------------------------------

    // Returns true if at least one sub-section with this name exists
    // (case-insensitive).
    bool has_section(const std::string& name) const {
        std::string lower = to_lower(name);
        for (const auto& kv : sections_)
            if (to_lower(kv.first) == lower) return true;
        return false;
    }

    // Returns a const reference to the first sub-section with this name
    // (case-insensitive). Throws if not found.
    const Config& section(const std::string& name) const {
        std::string lower = to_lower(name);
        for (const auto& kv : sections_)
            if (to_lower(kv.first) == lower) return kv.second;
        throw std::runtime_error("Section not found: " + name);
    }

    // Ordered list of all sub-sections (name, Config) pairs.
    // Preserves declaration order; may contain duplicate names.
    const std::vector<std::pair<std::string, Config>>& all_sections() const {
        return sections_;
    }

private:
    std::map<std::string, std::string>           values_;
    std::vector<std::pair<std::string, Config>>  sections_;  // ordered, allows duplicates

    static std::string to_lower(const std::string& s) {
        std::string r = s;
        for (auto& c : r) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
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
    size_t parse_block(const std::string& src, size_t pos) {
        std::string clean = remove_comments(src);
        return parse_block_inner(clean, pos);
    }

    size_t parse_block_inner(const std::string& src, size_t pos) {
        size_t n = src.size();

        auto skip_ws = [&](size_t p) {
            while (p < n && std::isspace(static_cast<unsigned char>(src[p]))) ++p;
            return p;
        };

        auto read_ident = [&](size_t p) -> std::pair<std::string, size_t> {
            std::string id;
            while (p < n && (std::isalnum(static_cast<unsigned char>(src[p]))
                             || src[p] == '_' || src[p] == '/' || src[p] == '.'))
                id += src[p++];
            return {id, p};
        };

        pos = skip_ws(pos);

        // If we're at top level, look for optional "Simulation {" wrapper
        if (pos == 0) {
            auto [name, after_name] = read_ident(pos);
            if (to_lower(name) == "simulation") {
                pos = skip_ws(after_name);
                if (pos < n && src[pos] == '{') {
                    ++pos;
                    pos = parse_inner(src, pos);
                    return pos;
                }
            }
            pos = 0;
        }

        return parse_inner(src, pos);
    }

    // Parse key=value pairs and Section { } blocks until '}' or end.
    size_t parse_inner(const std::string& src, size_t pos) {
        size_t n = src.size();

        auto skip_ws = [&](size_t p) {
            while (p < n && std::isspace(static_cast<unsigned char>(src[p]))) ++p;
            return p;
        };

        auto read_ident = [&](size_t p) -> std::pair<std::string, size_t> {
            std::string id;
            while (p < n && (std::isalnum(static_cast<unsigned char>(src[p]))
                             || src[p] == '_' || src[p] == '/' || src[p] == '.'))
                id += src[p++];
            return {id, p};
        };

        while (pos < n) {
            pos = skip_ws(pos);
            if (pos >= n) break;

            if (src[pos] == '}') {
                ++pos;
                pos = skip_ws(pos);
                if (pos < n && src[pos] == ';') ++pos;
                break;
            }

            auto [ident, after_ident] = read_ident(pos);
            if (ident.empty()) {
                ++pos;
                continue;
            }
            pos = skip_ws(after_ident);

            if (pos < n && src[pos] == '{') {
                ++pos;
                Config sub;
                pos = sub.parse_inner(src, pos);
                // push_back preserves order and allows duplicate names
                sections_.push_back({ident, std::move(sub)});
            } else if (pos < n && src[pos] == '=') {
                ++pos;
                pos = skip_ws(pos);
                std::string value;
                if (pos < n && src[pos] == '"') {
                    ++pos;
                    while (pos < n && src[pos] != '"') value += src[pos++];
                    if (pos < n) ++pos;
                } else {
                    while (pos < n && src[pos] != ';' && src[pos] != '\n')
                        value += src[pos++];
                    value = strip(value);
                }
                values_[ident] = value;
                pos = skip_ws(pos);
                if (pos < n && src[pos] == ';') ++pos;
            } else {
                while (pos < n && src[pos] != ';' && src[pos] != '{'
                                && src[pos] != '}') ++pos;
                if (pos < n && src[pos] == ';') ++pos;
            }
        }
        return pos;
    }
};

// ---------------------------------------------------------------------------
// get_or<T> specialisations
// ---------------------------------------------------------------------------

template<> inline double
Config::get_or<double>(const std::string& key, double dv) const {
    return get_double(key, dv);
}

template<> inline int
Config::get_or<int>(const std::string& key, int dv) const {
    return get_int(key, dv);
}

template<> inline std::string
Config::get_or<std::string>(const std::string& key, std::string dv) const {
    return get_string(key, dv);
}

template<> inline bool
Config::get_or<bool>(const std::string& key, bool dv) const {
    return get_bool(key, dv);
}

} // namespace atomistica
