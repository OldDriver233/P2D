#pragma once
#include <string>
#include <string_view>
#include <format>
#include <iostream>

enum class TokenType {
    END,
    EOL,
    IDENTIFIER,
    LITERAL,
    SCIENTIFIC,
    LEFT_PAREN,
    RIGHT_PAREN,
    PLUS,
    MINUS,
    ASTERISK,
    SLASH,
    BACK_SLASH,
    CARET,
};

class Token {
public:
    std::string_view sv;
    TokenType type;

    Token(std::string_view sv, TokenType type): sv(sv), type(type) {}

    void format() {
        if(type == TokenType::EOL) {
            //std::cout<<std::format("[EOL] {}\n", static_cast<int>(type));
        }
        else if(type == TokenType::END) {
            //std::cout<<std::format("[END] {}\n", static_cast<int>(type));
        }
        else {
            //std::cout<<std::format("{} {}\n", sv, static_cast<int>(type));
        }
    }
};