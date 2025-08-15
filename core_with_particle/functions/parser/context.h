#ifndef CONTEXT_H
#define CONTEXT_H

#include <map>
#include <llvm/Transforms/Scalar/LoopPassManager.h>

#undef I
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Module.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Passes/StandardInstrumentations.h"
#include "llvm/Transforms/InstCombine/InstCombine.h"
#include "llvm/Transforms/Scalar/GVN.h"
#include "llvm/Transforms/Scalar/Reassociate.h"
#include "llvm/Transforms/Scalar/SimplifyCFG.h"
#include "llvm/Transforms/IPO/FunctionAttrs.h"
#include "llvm/Transforms/Scalar/EarlyCSE.h"
#include "llvm/Transforms/Scalar/LoopInstSimplify.h"
#include "llvm/Transforms/Scalar/SROA.h"
#define I _Complex_I
#include "jit.h"

using namespace llvm;
using namespace llvm::orc;

struct CompileContext {
public:
    std::unique_ptr<LLVMContext> CContext;
    std::unique_ptr<IRBuilder<>> CBuilder;
    std::unique_ptr<Module> CModule;
    Function* CFunction;
    std::vector<Function*> CFunctions;
    std::vector<Value*> CVars;
    std::vector<Value*> CTmp;
    std::map<std::string, Function*> CExtern;
    std::unique_ptr<FunctionPassManager> TheFPM;
    std::unique_ptr<ModulePassManager> TheMPM;
    std::unique_ptr<LoopAnalysisManager> TheLAM;
    std::unique_ptr<FunctionAnalysisManager> TheFAM;
    std::unique_ptr<CGSCCAnalysisManager> TheCGAM;
    std::unique_ptr<ModuleAnalysisManager> TheMAM;
    std::unique_ptr<PassInstrumentationCallbacks> ThePIC;
    std::unique_ptr<StandardInstrumentations> TheSI;


    CompileContext(std::shared_ptr<orc::KaleidoscopeJIT> jit) {
        ExitOnError ExitOnErr;
        CContext = std::make_unique<LLVMContext>();
        CModule = std::make_unique<Module>("extern", *CContext);
        CBuilder = std::make_unique<IRBuilder<>>(*CContext);
        CModule->setDataLayout(jit->getDataLayout());

        std::vector<std::string> names = {"exp", "sin", "cos", "tan", "sinh", "cosh", "tanh"};
        for (auto x: names) {
            std::vector<Type *> Exts(1, Type::getDoubleTy(*CContext));
            FunctionType *FTs =
              FunctionType::get(Type::getDoubleTy(*CContext), Exts, false);
            Function* F =
                Function::Create(FTs, Function::ExternalLinkage, x, CModule.get());
            //CExtern.emplace(x, F);
        }

        std::vector<Type *> Exts(2, Type::getDoubleTy(*CContext));
        FunctionType *FTs =
          FunctionType::get(Type::getDoubleTy(*CContext), Exts, false);
        Function* F =
            Function::Create(FTs, Function::ExternalLinkage, "pow", CModule.get());
        //CExtern.emplace("pow", F);

        ExitOnErr(jit->addModule(
          ThreadSafeModule(std::move(CModule), std::move(CContext))));

        // We intentially repeats the process in order to let the jit know the common math functions.
        CContext = std::make_unique<LLVMContext>();
        CModule = std::make_unique<Module>("func", *CContext);
        CBuilder = std::make_unique<IRBuilder<>>(*CContext);
        CModule->setDataLayout(jit->getDataLayout());

        TheFPM = std::make_unique<FunctionPassManager>();
        TheMPM = std::make_unique<ModulePassManager>();
        TheLAM = std::make_unique<LoopAnalysisManager>();
        TheFAM = std::make_unique<FunctionAnalysisManager>();
        TheCGAM = std::make_unique<CGSCCAnalysisManager>();
        TheMAM = std::make_unique<ModuleAnalysisManager>();
        ThePIC = std::make_unique<PassInstrumentationCallbacks>();
        TheSI = std::make_unique<StandardInstrumentations>(*CContext,
                                                           /*DebugLogging*/ true);
        TheSI->registerCallbacks(*ThePIC, TheMAM.get());
        TheFPM->addPass(SimplifyCFGPass());
        TheFPM->addPass(SROAPass(SROAOptions::ModifyCFG));
        TheFPM->addPass(EarlyCSEPass());
        TheFPM->addPass(InstCombinePass());
        TheFPM->addPass(ReassociatePass());
        TheFPM->addPass(GVNPass());
        TheFPM->addPass(SimplifyCFGPass());

        PassBuilder PB;
        PB.registerModuleAnalyses(*TheMAM);
        PB.registerFunctionAnalyses(*TheFAM);
        PB.crossRegisterProxies(*TheLAM, *TheFAM, *TheCGAM, *TheMAM);

        for (auto x: names) {
            std::vector<Type *> Exts(1, Type::getDoubleTy(*CContext));
            FunctionType *FTs =
              FunctionType::get(Type::getDoubleTy(*CContext), Exts, false);
            Function* F =
                Function::Create(FTs, Function::ExternalLinkage, x, CModule.get());
            CExtern.emplace(x, F);
        }

        std::vector<Type *> Exts2(2, Type::getDoubleTy(*CContext));
        FunctionType *FTs2 =
          FunctionType::get(Type::getDoubleTy(*CContext), Exts2, false);
        Function* F2 =
            Function::Create(FTs2, Function::ExternalLinkage, "pow", CModule.get());
        CExtern.emplace("pow", F2);
    }

    void CreateFunction(int num_args, std::string name) {
        std::vector<Type *> Doubles(num_args, Type::getDoubleTy(*CContext));
        FunctionType *FT =
          FunctionType::get(Type::getDoubleTy(*CContext), Doubles, false);
        CFunction =
            Function::Create(FT, Function::ExternalLinkage, name, CModule.get());
        CFunctions.push_back(CFunction);

        CVars.clear();
        CTmp.clear();
        for (auto &arg: CFunction->args()) {
            CVars.push_back(&arg);
            CTmp.push_back(&arg);
        }

        BasicBlock *BB = BasicBlock::Create(*CContext, "entry", CFunction);
        CBuilder->SetInsertPoint(BB);
    }
};

#endif //CONTEXT_H
