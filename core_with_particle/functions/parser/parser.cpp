#include "parser.h"
#include "node.h"
#include "token.h"
#include <memory>
#include <stdexcept>

void Parser::tokenize() {
    std::string_view main = content;
    TokenType status = TokenType::END;
    std::size_t begin = 0;
    std::size_t line_cnt = 1, char_cnt = 0;
    bool ignore_eol = false;

    for (std::size_t i = 0; i < main.size(); i++) {
        char_cnt++;
        if (main[i] == ' ') {
            if (status != TokenType::END) {
                tokens.emplace_back(main.substr(begin, i - begin), status);
                status = TokenType::END;
            }
        } else if (main[i] == '\\') {
            ignore_eol = true;
        } else if (main[i] == '\n') {
            line_cnt++;
            char_cnt = 0;
            if (status != TokenType::END && !ignore_eol) {
                tokens.emplace_back(main.substr(begin, i - begin), status);
                tokens.emplace_back(main.substr(i, 1), TokenType::EOL);
                status = TokenType::END;
            }
            ignore_eol = false;
        } else {
            if(ignore_eol) {
                throw std::runtime_error(std::format("Unexpected '{}' after '\\'", main[i]));
            }
            if (main[i] == '+') {
                if (status != TokenType::END) {
                    tokens.emplace_back(main.substr(begin, i - begin), status);
                    status = TokenType::END;
                }
                tokens.emplace_back(main.substr(i, 1), TokenType::PLUS);
            }
            if (main[i] == '-') {
                if (status != TokenType::END) {
                    tokens.emplace_back(main.substr(begin, i - begin), status);
                    status = TokenType::END;
                }
                tokens.emplace_back(main.substr(i, 1), TokenType::MINUS);
            }
            if (main[i] == '*') {
                if (status != TokenType::END) {
                    tokens.emplace_back(main.substr(begin, i - begin), status);
                    status = TokenType::END;
                }
                tokens.emplace_back(main.substr(i, 1), TokenType::ASTERISK);
            }
            if (main[i] == '/') {
                if (status != TokenType::END) {
                    tokens.emplace_back(main.substr(begin, i - begin), status);
                    status = TokenType::END;
                }
                tokens.emplace_back(main.substr(i, 1), TokenType::SLASH);
            }
            if (main[i] == '^') {
                if (status != TokenType::END) {
                    tokens.emplace_back(main.substr(begin, i - begin), status);
                    status = TokenType::END;
                }
                tokens.emplace_back(main.substr(i, 1), TokenType::CARET);
            }
            if (main[i] == '(') {
                if (status != TokenType::END) {
                    tokens.emplace_back(main.substr(begin, i - begin), status);
                    status = TokenType::END;
                }
                tokens.emplace_back(main.substr(i, 1), TokenType::LEFT_PAREN);
            }
            if (main[i] == ')') {
                if (status != TokenType::END) {
                    tokens.emplace_back(main.substr(begin, i - begin), status);
                    status = TokenType::END;
                }
                tokens.emplace_back(main.substr(i, 1), TokenType::RIGHT_PAREN);
            }
            if (main[i] == '.') {
                if (status != TokenType::LITERAL) {
                    throw std::runtime_error(std::format(
                        "Error at {}:{}: Unexpected '.'", line_cnt, char_cnt));
                }
            }
            if (std::isdigit(main[i])) {
                if (status == TokenType::END) {
                    status = TokenType::LITERAL;
                    begin = i;
                }
            }
            if (std::isalpha(main[i])) {
                if (status == TokenType::END) {
                    status = TokenType::IDENTIFIER;
                    begin = i;
                }
                if (status == TokenType::LITERAL) {
                    throw std::runtime_error(
                        std::format("Error at {}:{}: Invalid literal value",
                                    line_cnt, char_cnt));
                }
            }
        }
    }
    if (status != TokenType::END) {
        tokens.emplace_back(main.substr(begin, main.size() - begin), status);
    }
    tokens.emplace_back("", TokenType::END);
}

Token Parser::peek() { return this->tokens[current]; }

Token Parser::peek_prev() { return this->tokens[current - 1]; }

Token Parser::advance() { return this->tokens[current++]; }

bool Parser::check(TokenType type) { return this->peek().type == type; }

bool Parser::match(TokenType type) {
    if (check(type)) {
        this->advance();
        return true;
    }
    return false;
}

Token Parser::consume(TokenType type, const std::string &msg) {
    if (check(type))
        return advance();
    throw std::runtime_error(msg);
}

void Parser::parse() { this->initial_node = this->term(); }

std::unique_ptr<Node> Parser::expr() { return this->term(); }

std::unique_ptr<Node> Parser::term() {
    auto expr = this->factor();
    while (match(TokenType::PLUS) || match(TokenType::MINUS)) {
        Token op = peek_prev();
        auto right = this->factor();
        if (op.type == TokenType::PLUS) {
            expr = std::make_unique<AddNode>(std::move(expr), std::move(right));
        } else {
            expr = std::make_unique<SubNode>(std::move(expr), std::move(right));
        }
    }
    return expr;
}

std::unique_ptr<Node> Parser::factor() {
    auto expr = this->unary();
    while (match(TokenType::ASTERISK) || match(TokenType::SLASH)) {
        Token op = peek_prev();
        auto right = this->unary();
        if (op.type == TokenType::ASTERISK) {
            expr =
                std::make_unique<MultNode>(std::move(expr), std::move(right));
        } else {
            expr = std::make_unique<DivNode>(std::move(expr), std::move(right));
        }
    }
    return expr;
}

std::unique_ptr<Node> Parser::unary() {
    if (match(TokenType::MINUS)) {
        auto expr = this->unary();
        return std::make_unique<NegateNode>(std::move(expr));
    }

    return this->call();
}

std::unique_ptr<Node> Parser::call() {
    auto expr = this->primary();

    if (check(TokenType::LEFT_PAREN)) {
        Token t = peek_prev();
        consume(TokenType::LEFT_PAREN, "Unreachable situation");
        expr = this->fetch_arg(t);
    }

    return expr;
}

std::unique_ptr<Node> Parser::fetch_arg(const Token &func_name) {
    [[unlikely]]
    if (func_name.type != TokenType::IDENTIFIER) {
        throw std::runtime_error("Expected identifier or operator before '('");
    } else {
        if (func_name.sv == "exp") {
            auto expr = this->expr();
            consume(TokenType::RIGHT_PAREN,
                    "Expected ')' after function calls");
            return std::make_unique<ExpNode>(std::move(expr));
        } else if (func_name.sv == "sin") {
            auto expr = this->expr();
            consume(TokenType::RIGHT_PAREN,
                    "Expected ')' after function calls");
            return std::make_unique<SinNode>(std::move(expr));
        } else if (func_name.sv == "cos") {
            auto expr = this->expr();
            consume(TokenType::RIGHT_PAREN,
                    "Expected ')' after function calls");
            return std::make_unique<CosNode>(std::move(expr));
        } else if (func_name.sv == "tan") {
            auto expr = this->expr();
            consume(TokenType::RIGHT_PAREN,
                    "Expected ')' after function calls");
            return std::make_unique<TanNode>(std::move(expr));
        } else if (func_name.sv == "sinh") {
            auto expr = this->expr();
            consume(TokenType::RIGHT_PAREN,
                    "Expected ')' after function calls");
            return std::make_unique<SinhNode>(std::move(expr));
        } else if (func_name.sv == "cosh") {
            auto expr = this->expr();
            consume(TokenType::RIGHT_PAREN,
                    "Expected ')' after function calls");
            return std::make_unique<CoshNode>(std::move(expr));
        } else if (func_name.sv == "tanh") {
            auto expr = this->expr();
            consume(TokenType::RIGHT_PAREN,
                    "Expected ')' after function calls");
            return std::make_unique<TanhNode>(std::move(expr));
        } else {
            throw std::runtime_error(
                std::format("Unknown builtin function name {}", func_name.sv));
        }
    }
}

std::unique_ptr<Node> Parser::primary() {
    if (match(TokenType::LITERAL)) {
        std::string tmp{this->peek_prev().sv};
        return std::make_unique<LiteralNode>(std::stod(tmp));
    }

    if (match(TokenType::LEFT_PAREN)) {
        auto expr = this->expr();
        consume(TokenType::RIGHT_PAREN,
                std::format("')' expected, got '{}'", peek().sv));
        return expr;
    }

    if (match(TokenType::IDENTIFIER)) {
        if (peek_prev().sv[0] == 'x') {
            std::string variable_id{peek_prev().sv.substr(1)};
            return std::make_unique<VariableNode>(std::stod(variable_id));
        }
        return std::make_unique<FuncNameNode>();
    }

    throw std::runtime_error(std::format("Unexpected '{}'", peek().sv));
};