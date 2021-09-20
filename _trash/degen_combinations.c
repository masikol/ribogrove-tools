// Pierre Lindenbaum is the author of this program: https://www.biostars.org/p/6219/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define RUN(base) copy[index]=base; recurs(seq,copy,index+1,len)

static void recurs(const char* seq, char* copy, int index, const int len) {
  if(index == len) {
    fwrite(copy, sizeof(char), len, stdout);
    fputc('\n', stdout);
  } else {
    switch(toupper(seq[index])) {
      case 'A': case 'T': case 'G': case 'C':  RUN(seq[index]); break;
      case 'N': RUN('A');RUN('T');RUN('G');RUN('C');break;
      case 'W': RUN('A');RUN('T');break;
      case 'S': RUN('G');RUN('C');break;
      case 'B': RUN('T');RUN('G');RUN('C');break;
      case 'D': RUN('A');RUN('T');RUN('G');break;
      case 'H': RUN('A');RUN('T');RUN('C');break;
      case 'V': RUN('A');RUN('G');RUN('C');break;
      case 'K': RUN('G');RUN('T');break;
      case 'M': RUN('A');RUN('C');break;
      case 'R': RUN('A');RUN('G');break;
      case 'Y': RUN('C');RUN('T');break;
      default: fprintf(stderr, "Bad base in %s (%c)\n", seq, seq[index]); exit(EXIT_FAILURE); break;
    }
  }
}

int main(int argc, char** argv) {
  char* seq;
  int len, i;

  if(argc != 2) {
    fprintf(stderr, "Usage : %s <dna>", argv[0]);
    return EXIT_FAILURE;
  }

  seq = argv[1];
  len = strlen(seq);

  char* copy=malloc((len+1)*sizeof(char));
  if(copy == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(EXIT_FAILURE);
  }
  copy[len] = '\0';

  recurs(seq, copy, 0, len);

  free(copy);

  return 0;
}
