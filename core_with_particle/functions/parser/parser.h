#pragma once
#include "token.h"
#include "node.h"
#include <cctype>
#include <format>
#include <fstream>
#include <ios>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class Parser {
public:
    std::string content;
    std::vector<Token> tokens;
    std::size_t current = 0;
    std::size_t vector_size;
    std::unique_ptr<Node> initial_node;

    Parser() = default;
    ~Parser() = default;

    Parser(std::string filename, std::size_t vector_size): vector_size(vector_size) {
        init(filename);
    }

    
    void init(std::string filename) {
        std::ifstream f(filename);
        std::stringstream ss;
        std::string line;
        while (getline(f, line)) {
            // std::cout<<line<<std::endl;
            ss << line << "\n";
        }

        f.close();
        content = ss.str();
        tokenize();
        parse();
    }
    void tokenize();
    void parse();
    std::unique_ptr<Node> expr();
    std::unique_ptr<Node> factor();
    std::unique_ptr<Node> power();
    std::unique_ptr<Node> unary();
    std::unique_ptr<Node> term();
    std::unique_ptr<Node> call();
    std::unique_ptr<Node> fetch_arg(const Token&);
    std::unique_ptr<Node> primary();
    Token peek();
    Token peek_prev();
    Token advance();
    bool match(TokenType type);
    Token consume(TokenType type, const std::string& msg);
    bool check(TokenType type);
};