#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>

void printUsage(int argc,char** argv){ printf("Fast UNIX string generator -- usage: %s <number of replicates> <length>\n",argv[0]); exit(-1); }

char *randstring(int length) {
    char const *string = "abcdefghijklmnopqrstuvwxyz0123456789";
    size_t stringLen = 26+10;
    char *randomString;

    randomString = (char*) malloc(sizeof(char)*(length+1));

    if (!randomString) {
        return (char*)0;
    }

    unsigned int key = 0;

    for (int n = 0;n < length;n++) {
        key = random() % stringLen;
        randomString[n] = string[key];
    }

    randomString[length] = '\0';
    return randomString;
}

int getUrandom(){
  int randomData = open("/dev/urandom", O_RDONLY);
  char myRandomData[25];
  size_t randomDataLen = 0;
  while (randomDataLen < sizeof myRandomData)
  {
      ssize_t result = read(randomData, myRandomData + randomDataLen, (sizeof myRandomData) - randomDataLen);
      if (result < 0)
      {
          printf("zero.\n");
          // error, unable to read /dev/random
      }
      randomDataLen += result;
  }
  close(randomData);
  // return the first digit we encounter as our random number seed
  for(int i=0;i<25;i++){
    if(isdigit(myRandomData[i])){
      return((int)myRandomData[i]);
    }
  }
  getUrandom(); /* this doesn't have a clean exit */
}

int main(int argc, char* argv[]){
  if(argc<2)
    printUsage(argc,argv);

  srandom(getUrandom());

  int replicates = (int) strtol(argv[1],NULL,10);
  int     length = (int) strtol(argv[2],NULL,10);

  for(int i=0;i<replicates;i++){
    char *s = (char*) malloc(sizeof(char) + 1);
    s=randstring(length);
    printf("%s\n",s);
  }
  exit(0);
}
