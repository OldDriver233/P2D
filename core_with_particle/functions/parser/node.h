#pragma once

#include <stdexcept>
#include <string>
#include <memory>
#include <iostream>
#include <format>

#include <eigen3/Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;

class Node {
public:
    virtual VectorXd eval(const Eigen::Ref<MatrixXd>&) = 0;
    virtual VectorXd eval_deriv(const Eigen::Ref<MatrixXd>&, int) = 0;
    virtual void show() = 0;
    ~Node() {}
};

class FuncNameNode: public Node {
public:
    FuncNameNode() {}
    ~FuncNameNode() {}
    FuncNameNode(const FuncNameNode& other) = default;

    VectorXd eval(const Eigen::Ref<MatrixXd>& _) override {
        (void) _;
        throw std::runtime_error("Internal error: this node should NOT exist in generated AST");
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& _, int __) override {
        (void) _;
        (void) __;
        throw std::runtime_error("Internal error: this node should NOT exist in generated AST");
    }

    void show() override {
        std::cout<<std::format("{{ FuncNameNode }}");
    }
};

class LiteralNode: public Node {
public:
    double value;
    LiteralNode() {}
    LiteralNode(double value) {
        this->value = value;
    }
    ~LiteralNode() {}
    LiteralNode(const LiteralNode& other) = default;

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(VectorXd::Ones(x.rows()) * value);
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int _) override {
        (void) _;
        return std::move(VectorXd::Zero(x.rows()));
    }

    void show() override {
        std::cout<<std::format("{{ LiteralNode value = {} }}", value);
    }
};

class VariableNode: public Node {
public:
    int variable_id;
    VariableNode() {}
    VariableNode(int id): variable_id(id) {}
    ~VariableNode() {}
    VariableNode(const VariableNode& other) = default;

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return x.col(variable_id);
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        if(variable_id == wrt) return std::move(VectorXd::Ones(x.rows()));
        else return std::move(VectorXd::Zero(x.rows()));
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(-val->eval(x));
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(-val->eval_deriv(x, wrt));
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(val->eval(x).array().exp());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(val->eval(x).array().exp() * val->eval_deriv(x, wrt).array());
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(val->eval(x).array().sin());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(val->eval(x).array().cos() * val->eval_deriv(x, wrt).array());
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(val->eval(x).array().cos());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(-val->eval(x).array().sin() * val->eval_deriv(x, wrt).array());
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(val->eval(x).array().tan());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        Eigen::ArrayXd cos = val->eval(x).array().cos().array();
        return std::move(val->eval_deriv(x, wrt).array() / (cos * cos));
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(val->eval(x).array().sinh());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(val->eval(x).array().cosh() * val->eval_deriv(x, wrt).array());
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(val->eval(x).array().cosh());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(val->eval(x).array().sinh() * val->eval_deriv(x, wrt).array());
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(val->eval(x).array().tanh());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        Eigen::ArrayXd cosh = val->eval(x).array().cosh().array();
        return std::move(val->eval_deriv(x, wrt).array() / (cosh * cosh));
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(left->eval(x).array() + right->eval(x).array());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(left->eval_deriv(x, wrt).array() + right->eval_deriv(x, wrt).array());
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(left->eval(x).array() - right->eval(x).array());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(left->eval_deriv(x, wrt).array() - right->eval_deriv(x, wrt).array());
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(left->eval(x).array() * right->eval(x).array());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        return std::move(left->eval_deriv(x, wrt).array() * right->eval(x).array() + left->eval(x).array() * right->eval_deriv(x, wrt).array());
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

    VectorXd eval(const Eigen::Ref<MatrixXd>& x) override {
        return std::move(left->eval(x).array() / right->eval(x).array());
    }

    VectorXd eval_deriv(const Eigen::Ref<MatrixXd>& x, int wrt) override {
        Eigen::ArrayXd rval = right->eval(x).array();
        return std::move((left->eval_deriv(x, wrt).array() * rval - left->eval(x).array() * right->eval_deriv(x, wrt).array()) / (rval * rval));
    }

    void show() override {
        std::cout<<std::format("{{ DivNode \nLeft = ");
        left->show();
        std::cout<<std::format("\nRight = ");
        right->show();
        std::cout<<std::format(" }}");
    }
};