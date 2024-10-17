#include "redis_connector.h"
#include <chrono>
#include <iostream>

redis_connector::redis_connector() {
    context = redisConnect("127.0.0.1", 6379);
    if (context == nullptr || context->err) {
        if (context) {
            std::cerr << "Error: " << context->errstr << std::endl;
            redisFree(context);
        } else {
            std::cerr << "Error: Can't allocate redis context" << std::endl;
        }
        exit(1);
    }
}

redis_connector::~redis_connector() {
    redisFree(context);
}

void redis_connector::rpush(const std::string &field, const std::string &value) {
    redisReply *reply = (redisReply *)redisCommand(context, "RPUSH %s %s", field.c_str(), value.c_str());
    if (reply == nullptr) {
        std::cerr << "Error: " << context->errstr << std::endl;
    }
    freeReplyObject(reply);
}

void redis_connector::set(const std::string &field, const std::string &value) {
    redisReply *reply = (redisReply *)redisCommand(context, "SET %s %s", field.c_str(), value.c_str());
    if (reply == nullptr) {
        std::cerr << "Error: " << context->errstr << std::endl;
    }
    freeReplyObject(reply);
}

void redis_connector::del(const std::string &field) {
    redisReply *reply = (redisReply *)redisCommand(context, "DEL %s", field.c_str());
    if (reply == nullptr) {
        std::cerr << "Error: " << context->errstr << std::endl;
    }
    freeReplyObject(reply);
}