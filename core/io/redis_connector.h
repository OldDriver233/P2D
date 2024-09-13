#ifndef FEM_REDIS_CONNECTOR
#define FEM_REDIS_CONNECTOR

#include <hiredis/hiredis.h>
#include <string>

class redis_connector {
public:
    redisContext *context;

    redis_connector();
    ~redis_connector();

    void rpush(const std::string& field, const std::string& value);
    void set(const std::string& field, const std::string& value);
    void del(const std::string& field);
};

#endif //FEM_REDIS_CONNECTOR