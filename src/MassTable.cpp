#include "MassTable.h"

void MassTable::Init(const char* file)
{
  int n,z,a;
  int counter=0;
  char sym[3];
  float massexcess;
  fFileName=file;
  //cout << "initializing mass table from file "<<file<<" ..."<<endl;

  std::ifstream input_file(file);
  if (!input_file){
    std::cout << "error opening mass table file"<< std::endl;
    exit(1);
  }

  while (input_file) {
    //    //cout << "before read"<<endl;
     input_file>>n;
    if (n >= 0) {
      input_file>>z>>a>>sym>>massexcess;
      fEntries[counter].N=n;
      fEntries[counter].Z=z;
      fEntries[counter].A=a;
      //      fEntries[counter].Symbol=sym;
      //      //cout <<"n,z,a,sym:" << n<<" "<<z<<" "<<a<<" "<<sym<<endl;
      //      strcpy(sym,fEntries[counter].Symbol);
      //      //cout <<"point 1"<<endl;
      sprintf(fEntries[counter].Symbol,"%s",sym);
      fEntries[counter].MassExcess=massexcess;
      //      //cout <<"point 2"<<endl;
      counter++;
      //      //cout << "end of loop"<<endl;
    }
  }
  //cout << "number of entries= "<<counter<<endl;
  fNumberOfEntries = counter;
  
  //cout << "mass table read with " << fNumberOfEntries <<" entries."<<endl ;
  input_file.close();
  //cout << "file "<<file<< " closed."<<endl;

  this->initstatus = true;

};

void MassTable::Print(void)
{
  for (int i=0; i<fNumberOfEntries; i++)
    {
      PrintEntry(i);
    }
};

void MassTable::PrintEntry(int i)
{
  //cout << "(N,Z,A)=("<<fEntries[i].N<<","<<fEntries[i].Z<<","<<
  // fEntries[i].A<<") Symbol = "<<fEntries[i].Symbol<<
	// 		   "Mass excess= "<<fEntries[i].MassExcess<<endl;
};

int MassTable::Lookup(int n, int z, int a)
{
  int counter=0;
  int returnvalue=0;
  //  //cout << "n,z,a= "<<n<<" "<<z<<" "<<a<<endl;
  if (a==0) a=n+z;
  for (counter=0;counter<fNumberOfEntries;counter++) {
    if (fEntries[counter].N==n && fEntries[counter].Z==z &&
	fEntries[counter].A==a) {
      returnvalue=1;
      return (returnvalue);
    }
  }
  return (returnvalue);
};

int MassTable::Lookup(char *symbol, int a)
{
  char test[3];
  sprintf(test,"%s",symbol);
  int counter=0;
  int returnvalue=0;
  //  //cout << "symbol= "<<symbol<<endl;
  for (counter=0;counter<fNumberOfEntries;counter++) {
    if (!strcmp(fEntries[counter].Symbol,test) && fEntries[counter].A==a){
      returnvalue=1;
      return (returnvalue);
    }
  }
  return (returnvalue);
};

float MassTable::GetMassExcess(const char *symbol, int a)
{
  char test[3];
  sprintf(test,"%s",symbol);
  int counter=0;
  float returnvalue=0;
  for (counter=0;counter<fNumberOfEntries;counter++) {
    if (!strcmp(fEntries[counter].Symbol,test) && fEntries[counter].A==a){
      returnvalue=fEntries[counter].MassExcess;
      return (returnvalue);
    }
  }
  return (returnvalue);
};

float MassTable::GetMassExcess(int n=0, int z=0, int a=0)
{
  int counter=0;
  float returnvalue=0;
  if (a==0) a=n+z;
  for (counter=0;counter<fNumberOfEntries;counter++) {
    if (fEntries[counter].N==n && fEntries[counter].Z==z &&
	fEntries[counter].A==a) {
      returnvalue=fEntries[counter].MassExcess;
      return (returnvalue);
    }
  }
  return (returnvalue);
};

const char* MassTable::GetSymbol(int z=0)
{
  int counter=0;
  const char* returnvalue="n ";
  for (counter=0;counter<fNumberOfEntries;counter++) {
    if (fEntries[counter].Z==z) {
      returnvalue=fEntries[counter].Symbol;
      return (returnvalue);
    }
  }
  return (returnvalue);
};

int MassTable::GetZ(char *symbol)
{
  char test[3];
  sprintf(test,"%s",symbol);
  int counter=0;
  int returnvalue=0;
  for (counter=0;counter<fNumberOfEntries;counter++) {
    if (!strcmp(fEntries[counter].Symbol,test)){
      returnvalue=fEntries[counter].Z;
      return (returnvalue);
    }
  }
  return (returnvalue);
};

//added by zmp:
bool MassTable::GetInitStatus(){//return true for initialized, false otherwise
  return this->initstatus;
};

float MassTable::GetMassMeV(const char* symbol, int a){//return the mass in MeV/c^2 instead of mass excess - simpler for the end user! (its me im the end user)
  //return (a*AMU + this->GetMassExcess(symbol,a)*0.001);

  if(a != 0 || (symbol != nullptr && strcmp(symbol, "None") == 0)){
    return 0.0f;
  }

  return (a*AMU + this->GetMassExcess(symbol,a)*0.001);

};
