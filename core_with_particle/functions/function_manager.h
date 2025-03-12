#ifndef FEM_FUNCTION_MANAGER_H
#define FEM_FUNCTION_MANAGER_H
#include "./parser/parser.h"
#include "../io/settings/settings.h"
#include "./parser/jit.h"
#include "./parser/context.h"
#include "llvm/Support/TargetSelect.h"
#include "llvm/Target/TargetMachine.h"
#include "llvm/Bitcode/BitcodeWriter.h"
#include "llvm/IR/Verifier.h"

class FunctionManager{
public:
    Parser uoc_anode;
    Parser uoc_cathode;
    Parser anode_entropy;
    Parser cathode_entropy;
    Parser kappa;
    Parser electrolyte_diffuse;
    double (*f_uoc_anode)(double) = nullptr;
    double (*f_d_uoc_anode)(double) = nullptr;
    double (*f_uoc_cathode)(double) = nullptr;
    double (*f_d_uoc_cathode)(double) = nullptr;
    double (*f_anode_entropy)(double) = nullptr;
    double (*f_cathode_entropy)(double) = nullptr;
    double (*f_kappa)(double, double) = nullptr;
    double (*f_d_kappa)(double, double) = nullptr;
    double (*f_diffuse_l)(double, double) = nullptr;
    double (*f_d_diffuse_l)(double, double) = nullptr;
    std::shared_ptr<KaleidoscopeJIT> TheJIT;

    FunctionManager() {
        ExitOnError ExitOnErr;
        InitializeNativeTarget();
        InitializeNativeTargetAsmPrinter();
        InitializeNativeTargetAsmParser();
        TheJIT = ExitOnErr(KaleidoscopeJIT::Create());
        CompileContext ctx(TheJIT);

        if(settings::use_customize_uoc) {
            uoc_anode.vector_size = 1;
            uoc_anode.init(settings::uoc_anode_path);
            uoc_cathode.vector_size = 1;
            uoc_cathode.init(settings::uoc_cathode_path);
            anode_entropy.vector_size = 1;
            anode_entropy.init(settings::anode_entropy_path);
            cathode_entropy.vector_size = 1;
            cathode_entropy.init(settings::cathode_entropy_path);

            ctx.CreateFunction(1, "u_anode");
            ctx.CBuilder->CreateRet(uoc_anode.initial_node->emit(ctx));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
            ctx.CreateFunction(1, "d_u_anode");
            ctx.CBuilder->CreateRet(uoc_anode.initial_node->emit_deriv(ctx, 0));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
            ctx.CreateFunction(1, "u_cathode");
            ctx.CBuilder->CreateRet(uoc_cathode.initial_node->emit(ctx));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
            ctx.CreateFunction(1, "d_u_cathode");
            ctx.CBuilder->CreateRet(uoc_cathode.initial_node->emit_deriv(ctx, 0));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
            ctx.CreateFunction(1, "dudt_anode");
            ctx.CBuilder->CreateRet(anode_entropy.initial_node->emit(ctx));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
            ctx.CreateFunction(1, "dudt_cathode");
            ctx.CBuilder->CreateRet(cathode_entropy.initial_node->emit(ctx));
            llvm::verifyFunction(*ctx.CFunction);
        }
        if (settings::use_customize_kappa) {
            kappa.vector_size = 2;
            kappa.init(settings::kappa_path);

            ctx.CreateFunction(2, "kappa");
            ctx.CBuilder->CreateRet(kappa.initial_node->emit(ctx));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
            ctx.CreateFunction(2, "d_kappa");
            ctx.CBuilder->CreateRet(kappa.initial_node->emit_deriv(ctx, 0));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
        }
        if (settings::use_customize_diffuse) {
            electrolyte_diffuse.vector_size = 2;
            electrolyte_diffuse.init(settings::diffuse_path);

            ctx.CreateFunction(2, "d_l");
            ctx.CBuilder->CreateRet(electrolyte_diffuse.initial_node->emit(ctx));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
            ctx.CreateFunction(2, "d_d_l");
            ctx.CBuilder->CreateRet(electrolyte_diffuse.initial_node->emit_deriv(ctx, 0));
            llvm::verifyFunction(*ctx.CFunction);
            ctx.TheFPM->run(*ctx.CFunction, *ctx.TheFAM);
        }

        auto RT = TheJIT->getMainJITDylib().createResourceTracker();
        auto TSM = ThreadSafeModule(std::move(ctx.CModule), std::move(ctx.CContext));
        ExitOnErr(TheJIT->addModule(std::move(TSM), RT));

        if(settings::use_customize_uoc) {
            auto ExprSymbol = ExitOnErr(TheJIT->lookup("u_anode"));
            f_uoc_anode = ExprSymbol.getAddress().toPtr<double (*)(double)>();
            ExprSymbol = ExitOnErr(TheJIT->lookup("d_u_anode"));
            f_d_uoc_anode = ExprSymbol.getAddress().toPtr<double (*)(double)>();
            ExprSymbol = ExitOnErr(TheJIT->lookup("u_cathode"));
            f_uoc_cathode = ExprSymbol.getAddress().toPtr<double (*)(double)>();
            ExprSymbol = ExitOnErr(TheJIT->lookup("d_u_cathode"));
            f_d_uoc_cathode = ExprSymbol.getAddress().toPtr<double (*)(double)>();
            ExprSymbol = ExitOnErr(TheJIT->lookup("dudt_anode"));
            f_anode_entropy = ExprSymbol.getAddress().toPtr<double (*)(double)>();
            ExprSymbol = ExitOnErr(TheJIT->lookup("dudt_cathode"));
            f_cathode_entropy = ExprSymbol.getAddress().toPtr<double (*)(double)>();
        }
        if (settings::use_customize_kappa) {
            auto ExprSymbol = ExitOnErr(TheJIT->lookup("kappa"));
            f_kappa = ExprSymbol.getAddress().toPtr<double (*)(double, double)>();
            ExprSymbol = ExitOnErr(TheJIT->lookup("d_kappa"));
            f_d_kappa = ExprSymbol.getAddress().toPtr<double (*)(double, double)>();
        }
        if (settings::use_customize_diffuse) {
            auto ExprSymbol = ExitOnErr(TheJIT->lookup("d_l"));
            f_diffuse_l = ExprSymbol.getAddress().toPtr<double (*)(double, double)>();
            ExprSymbol = ExitOnErr(TheJIT->lookup("d_d_l"));
            f_d_diffuse_l = ExprSymbol.getAddress().toPtr<double (*)(double, double)>();
        }
    }
};

#endif //FEM_FUNCTION_MANAGER_H