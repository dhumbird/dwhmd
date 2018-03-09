#include "sysfun.h"

time_t last=0;

string Hostname(){
  char hostname[256];
  gethostname(hostname, 256);
  string host=hostname;
  if (host.find(".")!=string::npos)
    host=host.replace(host.find("."), 1000,"");
  return host;
}

string DateTime(){
  time_t t1=time(0);
  string datetime=asctime(localtime(&t1));
  datetime=datetime.replace(datetime.find("\n"),1,"");
  return datetime;
}

string EstTimeFinish(int n){
  if (last==0){
    last=time(0);
    return "";
  }
  else{
    time_t newtime=time(0);
    newtime+=(int)difftime(newtime, last)*n;
    string donetime=asctime(localtime(&newtime));
    donetime=donetime.replace(donetime.find("\n"),1,"");
    last=newtime;
    return donetime;
  }
}

bool FileExists(string file){
  ifstream ifs(file.c_str(), ios::in);
  if (ifs) return 1;
  else return 0;
}

int system(char * s){
  int pid;
  
  if((pid=fork())==0){
    execlp("sh", "sh", "-c", s, (char *) 0);
    return pid;
  }
}
