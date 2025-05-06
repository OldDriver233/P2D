#pragma once

#include <stdexcept>
#include <string>
#include <memory>
#include <iostream>

#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Type.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/IR/Constants.h"
#include "context.h"

struct CompileContext;
using namespace llvm;

class Node {
public:
    virtual Value* emit(CompileContext &) = 0;
    virtual Value* emit_deriv(CompileContext&, int) = 0;
    virtual ~Node() {}
};

class FuncNameNode: public Node {
public:
    FuncNameNode() {}
    ~FuncNameNode() {}
    FuncNameNode(const FuncNameNode& other) = default;

    Value* emit(CompileContext &_) override {
        (void) _;
        //throw std::runtime_error("Internal error: this node should NOT exist in generated AST");
    }

    Value* emit_deriv(CompileContext &_, int __) override {
        (void) _;
        (void) __;
        //throw std::runtime_error("Internal error: this node should NOT exist in generated AST");
    }

};

class LiteralNode: public Node {
public:
    double value;
    LiteralNode() {}
    LiteralNode(double value) {
        this->value = value;
    }
    ~LiteralNode() = default;
    LiteralNode(const LiteralNode& other) = default;

    Value* emit(CompileContext &ctx) override {
        return ConstantFP::get(*ctx.CContext, APFloat(value));
    }

    Value* emit_deriv(CompileContext &ctx, int _) override {
        (void) _;
        return ConstantFP::get(*ctx.CContext, APFloat(0.0));
    }

};

class VariableNode: public Node {
public:
    int variable_id;
    VariableNode() {}
    VariableNode(int id): variable_id(id) {}
    ~VariableNode() {}
    VariableNode(const VariableNode& other) = default;

    Value* emit(CompileContext &ctx) override {
        return ctx.CTmp[variable_id];
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        if(variable_id == wrt) return ConstantFP::get(*ctx.CContext, APFloat(1.0));
        else return ConstantFP::get(*ctx.CContext, APFloat(0.0));
    }

};

class NegateNode: public Node {
public:
    std::unique_ptr<Node> val;
    NegateNode() {}
    NegateNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~NegateNode() {}
    NegateNode(const NegateNode& other) = delete;

    Value* emit(CompileContext &ctx) override {
        return ctx.CBuilder->CreateFNeg(val->emit(ctx));
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        return ctx.CBuilder->CreateFNeg(val->emit(ctx));
    }

};

class ExpNode: public Node {
public:
    std::unique_ptr<Node> val;
    ExpNode() {}
    ExpNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~ExpNode() {}
    ExpNode(const ExpNode& other) = delete;

    Value* emit(CompileContext &ctx) override {
        auto fn = ctx.CExtern.find("exp")->second;
        std::vector<Value*> V;
        V.push_back(val->emit(ctx));
        return ctx.CBuilder->CreateCall(fn, V);
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        auto fn = ctx.CExtern.find("exp")->second;
        std::vector<Value*> V, dV;
        V.push_back(val->emit(ctx));
        dV.push_back(val->emit_deriv(ctx, wrt));
        return ctx.CBuilder->CreateFMul(dV[0], ctx.CBuilder->CreateCall(fn, V));
    }

};

class SinNode: public Node {
public:
    std::unique_ptr<Node> val;
    SinNode() {}
    SinNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~SinNode() {}
    SinNode(const SinNode& other) = delete;

    Value* emit(CompileContext &ctx) override {
        auto fn = ctx.CExtern.find("sin")->second;
        std::vector<Value*> V;
        V.push_back(val->emit(ctx));
        return ctx.CBuilder->CreateCall(fn, V);
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        auto fn = ctx.CExtern.find("cos")->second;
        std::vector<Value*> V, dV;
        V.push_back(val->emit(ctx));
        dV.push_back(val->emit_deriv(ctx, wrt));
        return ctx.CBuilder->CreateFMul(dV[0], ctx.CBuilder->CreateCall(fn, V));
    }

};

class CosNode: public Node {
public:
    std::unique_ptr<Node> val;
    CosNode() {}
    CosNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~CosNode() {}
    CosNode(const CosNode& other) = delete;

    Value* emit(CompileContext &ctx) override {
        auto fn = ctx.CExtern.find("cos")->second;
        std::vector<Value*> V;
        V.push_back(val->emit(ctx));
        return ctx.CBuilder->CreateCall(fn, V);
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        auto fn = ctx.CExtern.find("sin")->second;
        std::vector<Value*> V, dV;
        V.push_back(val->emit(ctx));
        dV.push_back(val->emit_deriv(ctx, wrt));
        return ctx.CBuilder->CreateNeg(ctx.CBuilder->CreateFMul(dV[0], ctx.CBuilder->CreateCall(fn, V)));
    }

};

class TanNode: public Node {
public:
    std::unique_ptr<Node> val;
    TanNode() {}
    TanNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~TanNode() {}
    TanNode(const TanNode& other) = delete;

    Value* emit(CompileContext &ctx) override {
        auto fn = ctx.CExtern.find("tan")->second;
        std::vector<Value*> V;
        V.push_back(val->emit(ctx));
        return ctx.CBuilder->CreateCall(fn, V);
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        auto fn = ctx.CExtern.find("cos")->second;
        std::vector<Value*> V, dV;
        V.push_back(val->emit(ctx));
        dV.push_back(val->emit_deriv(ctx, wrt));
        return ctx.CBuilder->CreateFDiv(dV[0], ctx.CBuilder->CreateFMul(ctx.CBuilder->CreateCall(fn, V), ctx.CBuilder->CreateCall(fn, V)));
    }

};

class SinhNode: public Node {
public:
    std::unique_ptr<Node> val;
    SinhNode() {}
    SinhNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~SinhNode() {}
    SinhNode(const SinhNode& other) = delete;

    Value* emit(CompileContext &ctx) override {
        auto fn = ctx.CExtern.find("sinh")->second;
        std::vector<Value*> V;
        V.push_back(val->emit(ctx));
        return ctx.CBuilder->CreateCall(fn, V);
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        auto fn = ctx.CExtern.find("cosh")->second;
        std::vector<Value*> V, dV;
        V.push_back(val->emit(ctx));
        dV.push_back(val->emit_deriv(ctx, wrt));
        return ctx.CBuilder->CreateFMul(dV[0], ctx.CBuilder->CreateCall(fn, V));
    }

};

class CoshNode: public Node {
public:
    std::unique_ptr<Node> val;
    CoshNode() {}
    CoshNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~CoshNode() {}
    CoshNode(const CoshNode& other) = delete;

    Value* emit(CompileContext &ctx) override {
        auto fn = ctx.CExtern.find("cosh")->second;
        std::vector<Value*> V;
        V.push_back(val->emit(ctx));
        return ctx.CBuilder->CreateCall(fn, V);
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        auto fn = ctx.CExtern.find("sinh")->second;
        std::vector<Value*> V, dV;
        V.push_back(val->emit(ctx));
        dV.push_back(val->emit_deriv(ctx, wrt));
        return ctx.CBuilder->CreateFMul(dV[0], ctx.CBuilder->CreateCall(fn, V));
    }

};

class TanhNode: public Node {
public:
    std::unique_ptr<Node> val;
    TanhNode() {}
    TanhNode(std::unique_ptr<Node> val): val(std::move(val)) {}
    ~TanhNode() {}
    TanhNode(const TanhNode& other) = delete;

    Value* emit(CompileContext &ctx) override {
        auto fn = ctx.CExtern.find("tanh")->second;
        std::vector<Value*> V;
        V.push_back(val->emit(ctx));
        return ctx.CBuilder->CreateCall(fn, V);
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        auto fn = ctx.CExtern.find("cosh")->second;
        std::vector<Value*> V, dV;
        V.push_back(val->emit(ctx));
        dV.push_back(val->emit_deriv(ctx, wrt));
        return ctx.CBuilder->CreateFDiv(dV[0], ctx.CBuilder->CreateFMul(ctx.CBuilder->CreateCall(fn, V), ctx.CBuilder->CreateCall(fn, V)));
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

    Value* emit(CompileContext &ctx) override {
        return ctx.CBuilder->CreateFAdd(left->emit(ctx), right->emit(ctx));
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        return ctx.CBuilder->CreateFAdd(left->emit_deriv(ctx, wrt), right->emit_deriv(ctx, wrt));
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

    Value* emit(CompileContext &ctx) override {
        return ctx.CBuilder->CreateFSub(left->emit(ctx), right->emit(ctx));
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        return ctx.CBuilder->CreateFSub(left->emit_deriv(ctx, wrt), right->emit_deriv(ctx, wrt));
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
    Value* emit(CompileContext &ctx) override {
        return ctx.CBuilder->CreateFMul(left->emit(ctx), right->emit(ctx));
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        return ctx.CBuilder->CreateFAdd(
            ctx.CBuilder->CreateFMul(left->emit(ctx), right->emit_deriv(ctx, wrt)),
            ctx.CBuilder->CreateFMul(left->emit_deriv(ctx, wrt), right->emit(ctx))
            );
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

    Value* emit(CompileContext &ctx) override {
        return ctx.CBuilder->CreateFDiv(left->emit(ctx), right->emit(ctx));
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        return ctx.CBuilder->CreateFDiv(
                ctx.CBuilder->CreateFSub(
                    ctx.CBuilder->CreateFMul(left->emit_deriv(ctx, wrt), right->emit(ctx)),
                    ctx.CBuilder->CreateFMul(left->emit(ctx), right->emit_deriv(ctx, wrt))
                    ),
                ctx.CBuilder->CreateFMul(
                    right->emit(ctx),
                    right->emit(ctx)
                    )
            );
    }

};

class PowNode: public Node {
public:
    std::unique_ptr<Node> left, right;
    PowNode() {}
    PowNode(std::unique_ptr<Node> left, std::unique_ptr<Node> right): left(std::move(left)), right(std::move(right)) {}
    ~PowNode() {}
    PowNode(const DivNode& other) = delete;
    PowNode(PowNode&& other) = default;

    Value* emit(CompileContext &ctx) override {
        auto fn = ctx.CModule->getFunction("pow");
        std::vector<Value*> V;
        V.push_back(left->emit(ctx));
        V.push_back(right->emit(ctx));
        return ctx.CBuilder->CreateCall(fn, V);
    }

    Value* emit_deriv(CompileContext &ctx, int wrt) override {
        auto fn = ctx.CModule->getFunction("pow");
        std::vector<Value*> V, dV;
        V.push_back(left->emit(ctx));
        V.push_back(ctx.CBuilder->CreateFSub(right->emit(ctx), ConstantFP::get(*ctx.CContext, APFloat(1.0))));
        dV.push_back(left->emit_deriv(ctx, wrt));
        return ctx.CBuilder->CreateFMul(
            dV[0],
            ctx.CBuilder->CreateFMul(right->emit(ctx), ctx.CBuilder->CreateCall(fn, V)));
    }

};