// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskMyTask.h"
#include "TInterpreter.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "TMacro.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliMultSelectionTask.h"
#include "TSystem.h"
#include "TChain.h"
#include "AliAnalysisAlien.h"
#include "AliPhysicsSelectionTask.h"

//==================================================================
//注意！！！！！
//进入ali环境后，aliroot runAnalysis.C 执行此脚本提交任务
//在local运行时，设置 Bool_t local = kTRUE ，Bool gridTest = kTRUE
//在Grid上运行时，设置 Bool_t local = kFALSE，Bool gridTest = kFALSE
//有一个ROOT5没有的模式称之为GridTest，可以不提交grid任务而进行远程测试，设置 Bool_t local = kFALSE，Bool gridTest = kTRUE


//当Grid上任务完成，不用改变 Bool_t local = kFALSE，Bool gridTest = kFALSE 的参数设置，
//把213行的 alienHandler->SetRunMode("full") 改为 alienHandler->SetRunMode("terminate")
//再在终端输入 aliroot runAnalysis.C 就会提交在run内的merge任务

//run内merge完之后，将193行的 alienHandler->SetMergeViaJDL(kTRUE) 改为 alienHandler->SetMergeViaJDL(kFALSE)
//其他部分均不要改变
//再在终端输入 aliroot runAnalysis.C 就会提交在run之间的merge任务，直接会将结果下载到本地，不需要去Monalisa下载
//==================================================================


void runAnalysis(int runMode = 0)
{
    // Dataset
    TString dataset = "LHC15o";//设置dataset名称
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
   Bool_t local;
Bool_t gridTest;
// mode 0 = local test
// mode 1 = grid test
// mode 2 = grid run merge
// mode 3 = run 之间 merge
if(runMode==0){
local = kTRUE;
gridTest = kTRUE;
}
if(runMode==1||runMode==2||runMode==3){
local = kFALSE;
gridTest = kFALSE;
}
	 //Bool_t local = kFALSE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    //Bool_t gridTest = kFALSE;
    

/*
	Bool_t local = kTRUE;
	Bool_t gridTest = kTRUE;
*/
    // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
     
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
    AliAnalysisTaskPIDResponse* PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse*>(PIDadd.Exec());
   
    if (dataset.EqualTo("LHC15o") || dataset.EqualTo("LHC18q") || dataset.EqualTo("LHC18r")) {
      TMacro physicsSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));
      AliPhysicsSelectionTask* physicsSelectionTask = reinterpret_cast<AliPhysicsSelectionTask*>(physicsSelection.Exec());
      TMacro multSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
      AliMultSelectionTask* multSelectionTask = reinterpret_cast<AliMultSelectionTask*>(multSelection.Exec());
    }
    
    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskMyTask.cxx++g");
    AliAnalysisTaskMyTask *task = reinterpret_cast<AliAnalysisTaskMyTask*>(gInterpreter->ExecuteMacro("AddMyTask.C"));
#else
    gROOT->LoadMacro("AliAnalysisTaskMyTask.cxx++g");
    gROOT->LoadMacro("AddMyTask.C");
    AliAnalysisTaskMyTask *task = AddMyTask();
#endif


    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)
/*
        if(dataset.EqualTo("LHC10h")) chain->Add("/Users/wangchunzheng/alice/data/2010/LHC10h/000139510/ESDs/pass2/AOD160/0247/AliAOD.root");//这里要改成本地Data文件路径
        if(dataset.EqualTo("LHC15o")) chain->Add("/Users/wangchunzheng/alice/data/2015/LHC15o/000245151/pass2/AOD252/0008/AliAOD.root");
        if(dataset.EqualTo("LHC18q")) chain->Add("/Users/wangchunzheng/alice/data/2018/LHC18q/000295588/pass3/AOD252/AOD/001/AliAOD.root");
        if(dataset.EqualTo("LHC18r")) chain->Add("/Users/wangchunzheng/alice/data/2018/LHC18r/000296691/pass3/AOD252/0001/AliAOD.root");
        // start the analysis locally, reading the events from the tchain
  
*/
//chain->Add("/home/df/work/2021Test/Pass2_AliAOD.root");
 chain->Add("/afs/cern.ch/user/d/dowang/public/AliAOD.root"); 
    mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskMyTask.cxx AliAnalysisTaskMyTask.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskMyTask.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20220420_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        if (dataset.EqualTo("LHC10h")) {
          //10h
          // select the input data
          alienHandler->SetGridDataDir("/alice/data/2010/LHC10h");
          alienHandler->SetDataPattern("ESDs/pass2/AOD160/*/AliAOD.root");
          // MC has no prefix, data has prefix 000
          alienHandler->SetRunPrefix("000");
          // runnumber
          // alienHandler->AddRunNumber(139510);
          // alienHandler->AddRunNumber(139507);
          // alienHandler->AddRunNumber(139505);
          // alienHandler->AddRunNumber(139503);
          // alienHandler->AddRunNumber(139465);
          // alienHandler->AddRunNumber(139438);
          // alienHandler->AddRunNumber(139437);
          // alienHandler->AddRunNumber(139360);
          // alienHandler->AddRunNumber(139329);
          // alienHandler->AddRunNumber(139328);
          alienHandler->AddRunNumber(139314);
          alienHandler->AddRunNumber(139310);
        }
        if (dataset.EqualTo("LHC15o")) {
          // select the input data
          alienHandler->SetGridDataDir("/alice/data/2015/LHC15o");
          alienHandler->SetDataPattern("pass2/AOD252/AOD/*/AliAOD.root");
          // MC has no prefix, data has prefix 000
          alienHandler->SetRunPrefix("000");
          alienHandler->AddRunNumber(244917);
          // alienHandler->AddRunNumber(244918);
          // alienHandler->AddRunNumber(244975);
          // alienHandler->AddRunNumber(244980);
          // alienHandler->AddRunNumber(244982);
          // alienHandler->AddRunNumber(244983);
          // alienHandler->AddRunNumber(245061);
          // alienHandler->AddRunNumber(245064);
        }
        if (dataset.EqualTo("LHC18q")) {
          // select the input data
          alienHandler->SetGridDataDir("/alice/data/2018/LHC18q");
          alienHandler->SetDataPattern("pass3/AOD252/AOD/*/AliAOD.root");
          // MC has no prefix, data has prefix 000
          alienHandler->SetRunPrefix("000");
          alienHandler->AddRunNumber(296269);
          alienHandler->AddRunNumber(296270); 
          alienHandler->AddRunNumber(296547);
          alienHandler->AddRunNumber(296472);
          alienHandler->AddRunNumber(295581);
          alienHandler->AddRunNumber(295584);
          alienHandler->AddRunNumber(295585);
          alienHandler->AddRunNumber(295586);
        }
        if (dataset.EqualTo("LHC18r")) {
          // select the input data
          alienHandler->SetGridDataDir("/alice/data/2018/LHC18r");
          alienHandler->SetDataPattern("pass3/AOD252/AOD/*/AliAOD.root");
          // MC has no prefix, data has prefix 000
          alienHandler->SetRunPrefix("000");
          alienHandler->AddRunNumber(296690);
          alienHandler->AddRunNumber(296691);
          alienHandler->AddRunNumber(296693);
          alienHandler->AddRunNumber(296694);
          alienHandler->AddRunNumber(296752);
          alienHandler->AddRunNumber(296784);
          alienHandler->AddRunNumber(296785);
          alienHandler->AddRunNumber(296786);
          alienHandler->AddRunNumber(296787);
        }

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        if(runMode==1||runMode==2) alienHandler->SetMergeViaJDL(kTRUE);
        if(runMode==3) alienHandler->SetMergeViaJDL(kFALSE);

        // define the output folders
        if (dataset.EqualTo("LHC10h")) alienHandler->SetGridWorkingDir("DCACheck/LHC10h");
        if (dataset.EqualTo("LHC15o")) alienHandler->SetGridWorkingDir("DCACheck/LHC15o");
        if (dataset.EqualTo("LHC18q")) alienHandler->SetGridWorkingDir("DCACheck/LHC18q");
        if (dataset.EqualTo("LHC18r")) alienHandler->SetGridWorkingDir("DCACheck/LHC18r");
        alienHandler->SetGridOutputDir("output");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            if(runMode==1) alienHandler->SetRunMode("full");
            if(runMode==2||runMode==3) alienHandler->SetRunMode("terminate");
            mgr->StartAnalysis("grid");
        }
    }
}

