#ifndef MASSTABLE_H
#define MASSTABLE_H

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <string.h>

struct MassTableEntry
{
  int N;
  int Z;
  int A;
  char Symbol[3];
  float MassExcess;
};
  
class MassTable
{
private:
  const char* fFileName;
  MassTableEntry fEntries[4000];
  int fNumberOfEntries;
  //bool initstatus = false;
  const float AMU=931.494013;

public:
  //  TMassTable();
  void Init(const char* file);
  int Lookup(int n, int z, int a);
  int Lookup(char* symbol, int a);
  float GetMassExcess(int n, int z, int a);
  float GetMassExcess(const char* symbol, int a);
  const char* GetSymbol(int z);
  int GetZ(char* symbol);
  void Print(void);
  void PrintEntry(int);
  //added by zmp
  bool GetInitStatus();
  bool initstatus = false;
  float GetMassMeV(const char* symbol, int a);
};

#endif