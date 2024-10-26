#pragma once

#include <stdexcept>
#include <string>
#include <memory>
#include <iostream>
#include <format>

#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
using namespace autodiff;

class Node {
public:
    virtual real eval(const Eigen::Ref<VectorXreal>&) = 0;
    virtual void show() = 0;
    ~Node() {}
};

class FuncNameNode: public Node {
public:
    FuncNameNode() {}
    ~FuncNameNode() {}
    FuncNameNode(const FuncNameNode& other) = default;

    real eval(const Eigen::Ref<VectorXreal>& _) override {
        (void) _;
        throw std::runtime_error("Internal error: this node should NOT exist in generated AST");
    }

    void show() override {
        std::cout<<std::format("{{ FuncNameNode }}");
    }
};

class LiteralNode: public Node {
public:
    real value;
    LiteralNode() {}
    LiteralNode(double value) {
        this->value = value;
    }
    ~LiteralNode() {}
    LiteralNode(const LiteralNode& other) = default;

    real eval(const Eigen::Ref<VectorXreal>& _) override {
        (void) _;
        return value;
    }

    void show() override {
        std::cout<<std::format("{{ LiteralNode value = {} }}", value.val());
    }
};

class VariableNode: public Node {
public:
    int variable_id;
    VariableNode() {}
    VariableNode(int id): variable_id(id) {}
    ~VariableNode() {}
    VariableNode(const VariableNode& other) = default;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return x(variable_id);
    }

    void show() override {
        std::cout<<std::format("{{ VariableNode variable_id = {} }}", variable_id);
    }
};

class NegateNode: public Node {
public:
    std::unique_ptr<Node> val;
    NegateNode() {}
    NegateNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~NegateNode() {}
    NegateNode(const NegateNode& other) = delete;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return -val->eval(x);
    }

    void show() override {
        std::cout<<std::format("{{ NegateNode \n value = ");
        val->show();
        std::cout<<std::format("\n}}");
    }
};

class ExpNode: public Node {
public:
    std::unique_ptr<Node> val;
    ExpNode() {}
    ExpNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~ExpNode() {}
    ExpNode(const ExpNode& other) = delete;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return exp(val->eval(x));
    }

    void show() override {
        std::cout<<std::format("{{ ExpNode \n value = ");
        val->show();
        std::cout<<std::format("\n}}");
    }
};

class SinNode: public Node {
public:
    std::unique_ptr<Node> val;
    SinNode() {}
    SinNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~SinNode() {}
    SinNode(const SinNode& other) = delete;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return sin(val->eval(x));
    }

    void show() override {
        std::cout<<std::format("{{ SinNode \n value = ");
        val->show();
        std::cout<<std::format("\n}}");
    }
};

class CosNode: public Node {
public:
    std::unique_ptr<Node> val;
    CosNode() {}
    CosNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~CosNode() {}
    CosNode(const CosNode& other) = delete;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return cos(val->eval(x));
    }

    void show() override {
        std::cout<<std::format("{{ CosNode \n value = ");
        val->show();
        std::cout<<std::format("\n}}");
    }
};

class TanNode: public Node {
public:
    std::unique_ptr<Node> val;
    TanNode() {}
    TanNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~TanNode() {}
    TanNode(const TanNode& other) = delete;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return tan(val->eval(x));
    }

    void show() override {
        std::cout<<std::format("{{ TanNode \n value = ");
        val->show();
        std::cout<<std::format("\n}}");
    }
};

class SinhNode: public Node {
public:
    std::unique_ptr<Node> val;
    SinhNode() {}
    SinhNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~SinhNode() {}
    SinhNode(const SinhNode& other) = delete;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return sinh(val->eval(x));
    }

    void show() override {
        std::cout<<std::format("{{ SinhNode \n value = ");
        val->show();
        std::cout<<std::format("\n}}");
    }
};

class CoshNode: public Node {
public:
    std::unique_ptr<Node> val;
    CoshNode() {}
    CoshNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~CoshNode() {}
    CoshNode(const CoshNode& other) = delete;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return cosh(val->eval(x));
    }

    void show() override {
        std::cout<<std::format("{{ CoshNode \n value = ");
        val->show();
        std::cout<<std::format("\n}}");
    }
};

class TanhNode: public Node {
public:
    std::unique_ptr<Node> val;
    TanhNode() {}
    TanhNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~TanhNode() {}
    TanhNode(const TanhNode& other) = delete;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return tanh(val->eval(x));
    }

    void show() override {
        std::cout<<std::format("{{ TanhNode \n value = ");
        val->show();
        std::cout<<std::format("\n}}");
    }
};

class AddNode: public Node {
public:
    std::unique_ptr<Node> left, right;
    AddNode() {}
    AddNode(std::unique_ptr<Node> left, std::unique_ptr<Node> right): left(std::move(left)), right(std::move(right)) {}
    ~AddNode() {}
    AddNode(const AddNode& other) = delete;
    AddNode(AddNode&& other) = default;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return left->eval(x) + right->eval(x);
    }

    void show() override {
        std::cout<<std::format("{{ AddNode \nLeft = ");
        left->show();
        std::cout<<std::format("\nRight = ");
        right->show();
        std::cout<<std::format("\n}}");
    }
};

class SubNode: public Node {
public:
    std::unique_ptr<Node> left, right;
    SubNode() {}
    SubNode(std::unique_ptr<Node> left, std::unique_ptr<Node> right): left(std::move(left)), right(std::move(right)) {}
    ~SubNode() {}
    SubNode(const SubNode& other) = delete;
    SubNode(SubNode&& other) = default;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return left->eval(x) - right->eval(x);
    }

    void show() override {
        std::cout<<std::format("{{ SubNode \nLeft = ");
        left->show();
        std::cout<<std::format("\nRight = ");
        right->show();
        std::cout<<std::format(" }}");
    }
};

class MultNode: public Node {
public:
    std::unique_ptr<Node> left, right;
    MultNode() {}
    MultNode(std::unique_ptr<Node> left, std::unique_ptr<Node> right): left(std::move(left)), right(std::move(right)) {}
    ~MultNode() {}
    MultNode(const MultNode& other) = delete;
    MultNode(MultNode&& other) = default;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return left->eval(x) * right->eval(x);
    }

    void show() override {
        std::cout<<std::format("{{ MultNode \nLeft = ");
        left->show();
        std::cout<<std::format("\nRight = ");
        right->show();
        std::cout<<std::format(" }}");
    }
};

class DivNode: public Node {
public:
    std::unique_ptr<Node> left, right;
    DivNode() {}
    DivNode(std::unique_ptr<Node> left, std::unique_ptr<Node> right): left(std::move(left)), right(std::move(right)) {}
    ~DivNode() {}
    DivNode(const DivNode& other) = delete;
    DivNode(DivNode&& other) = default;

    real eval(const Eigen::Ref<VectorXreal>& x) override {
        return left->eval(x) / right->eval(x);
    }

    void show() override {
        std::cout<<std::format("{{ DivNode \nLeft = ");
        left->show();
        std::cout<<std::format("\nRight = ");
        right->show();
        std::cout<<std::format(" }}");
    }
};