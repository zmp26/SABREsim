#include <iostream>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>
#include <TH1.h>
#include <TCutG.h>
#include <TGraph.h>
#include <TTree.h>
#include <TCanvas.h>

void RecursiveCopyDir(TDirectory *source, TDirectory *target) {
	// first, copy objects stored in ROOT's memory list for this directory
	TList *memList = source->GetList();
	if (memList) {
		TIter nextMem(memList);
		TObject *obj = nullptr;
		while ((obj = nextMem())) {
			if (!obj) continue;

			// skip file handles to prevent infinite recursion
			if (obj->InheritsFrom(TFile::Class())) continue;

			target->cd();
			
			// TTrees bound to an existing input file must be cloned to the new file target
			if (obj->InheritsFrom(TTree::Class())) {
				TTree *oldTree = dynamic_cast<TTree*>(obj);
				if (oldTree) {
					TTree *newTree = oldTree->CloneTree(-1, "fast");
					newTree->Write();
				}
			} else {
				if(obj->InheritsFrom(TH1::Class())) {
					dynamic_cast<TH1*>(obj)->SetOption("HIST");
				}
				obj->Write(obj->GetName(), TObject::kOverwrite);
			}
		}
	}

	// next, copy the keys (handles existing subdirectories and un-instantiated keys)
	TList *keyList = source->GetListOfKeys();
	if (keyList) {
		TIter nextKey(keyList);
		TKey *key = nullptr;
		while ((key = dynamic_cast<TKey*>(nextKey()))) {
			TObject *obj = key->ReadObj();
			if (!obj) continue;

			if (obj->InheritsFrom(TDirectory::Class())) {
				TDirectory *srcSubDir = dynamic_cast<TDirectory*>(obj);
				target->cd();
				TDirectory *tarSubDir = target->mkdir(srcSubDir->GetName());
				RecursiveCopyDir(srcSubDir, tarSubDir);
			} else if (obj->InheritsFrom(TTree::Class())) {
				//check if tree was already written via memory list to avoid duplicate entries
				target->cd();
				if (!target->Get(key->GetName())) {
					TTree *oldTree = dynamic_cast<TTree*>(obj);
					if (oldTree) {
						TTree *newTree = oldTree->CloneTree(-1, "fast");
						newTree->Write();
					}
				}
			} else {
				target->cd();
				obj->Write(key->GetName(), TObject::kOverwrite);
			}
		}
	}
}

void SaveSession(const char* outputFilename = "shared_analysis_results.root") {
	//save current working directory reference
	TDirectory *currentDir = gDirectory;
	if (!currentDir) {
		std::cerr << "Error: No active ROOT directory found." << std::endl;
		return;
	}

	std::cout << "Creating package ROOT file: " << outputFilename << "..." << std::endl;

	TFile *outFile = TFile::Open(outputFilename, "RECREATE");
	if (!outFile || outFile->IsZombie()) {
		std::cerr << "Error: Could not open output file " << outputFilename << " for writing." << std::endl;
		return;
	}

	//recursively export all objects and subdirectories into the destination file
	RecursiveCopyDir(currentDir, outFile);

	outFile->cd();
	auto *specials = gROOT->GetListOfSpecials();
	if(specials){
		TIter nextSpecial(specials);
		TObject *obj = nullptr;
		int cutCount = 0;
		while((obj = nextSpecial())){
			if(obj && obj->InheritsFrom(TCutG::Class())){
				obj->Write(obj->GetName(), TObject::kOverwrite);
				cutCount++;
			}
		}
		if(cutCount > 0) std::cout << "Saved " << cutCount << " TCutG object(s) from gROOT specials." << std::endl;
	}

	outFile->Close();
	delete outFile;

	currentDir->cd();
	std::cout << "Session successfully exported to '" << outputFilename << "'!" << std::endl;
}