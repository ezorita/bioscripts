#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
   char colspace [4][4] = {
      {0,1,2,3},
      {1,0,3,2},
      {2,3,0,1},
      {3,2,1,0}
   };

   char to_base[4] = {'A','C','G','T'};

   // Check params.
   if (argc != 2) {
      fprintf(stderr, "usage: %s <FASTA/FASTQ color sequence file>\n",argv[0]);
      exit(1);
   }

   // Input file format.
   FILE * fin = fopen(argv[1],"r");
   if (fin == NULL) {
      fprintf(stderr, "error: opening input file 'fopen()'.\n");
      exit(1);
   }
   int call_zcat = 0;
   unsigned char c = fgetc(fin);
   if (c == 0x1f) {
      c = fgetc(fin);
      if (c == 0x8b) {
         call_zcat = 1;
      }
      ungetc(c,fin);
      ungetc(0x1f,fin);
   } else ungetc(c,fin);

   // Fork and decompress using 'zcat'.
   int pipefd[2];
   if (call_zcat) {
      fclose(fin);
      if (pipe(pipefd) < 0) {
         fprintf(stderr, "error: 'pipe()'\n");
         exit(1);
      }
      int pid = fork();
      if (pid == -1) {
         fprintf(stderr, "error calling 'zcat': 'fork()'\n");
         fprintf(stderr, "error: the input stream seems to be compressed (gzip).\n");
         fprintf(stderr, "uncompress it first or use:\n  zcat %s | %s /dev/stdin\n",argv[1], argv[0]);
         exit(1);
      } else if (pid == 0) {
         close(pipefd[0]);
         // Replace stdout by pipe write end.
         close(1);
         if (dup(pipefd[1]) < 0) {
            fprintf(stderr, "error: 'dup()'\n");
            close(pipefd[1]);
            exit(1);
         }
         // Execv 'zcat'.
         char * args[3] = {"zcat",argv[1],NULL};
         if (execv("/bin/zcat", args) < 0) {
            fprintf(stderr, "error calling 'zcat': 'execv()' (is zcat in $PATH?)\n");
            fprintf(stderr, "error: the input stream seems to be compressed (gzip).\n");
            fprintf(stderr, "uncompress it first or use:\n  zcat %s | %s /dev/stdin\n",argv[1], argv[0]);
            close(pipefd[1]);
            exit(1);
         }
         exit(0);
      }
      // Close pipe write end.
      close(pipefd[1]);
      fin = fdopen(pipefd[0], "r");
   }


   // Detect file format.
   c = fgetc(fin);
   int div, seqline, qline;
   if (c == '>') {
      div = 2; seqline = 1; qline = 1;
   } else if (c == '@') {
      div = 4; seqline = 1; qline = 3;
   } else if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
      div = 1; seqline = 0; qline = 0;
   } else {
      fprintf(stderr, "error: file format not recognized.\n");
      fprintf(stderr, "accepted formats: FASTA, FASTQ or RAW.\n");
      fclose(fin);
      exit(1);
   }
   ungetc(c,fin);

   // Read input file.
   size_t  bufsize = 100;
   ssize_t bytesrd;
   char *  line = malloc(bufsize);
   size_t  lineno = 0;

   while((bytesrd = getline(&line, &bufsize, fin)) > 0) {
      if (lineno % div == seqline) {
         char * bases = malloc(bytesrd-1);
         bases[bytesrd-2] = 0;
         // Read reference base.
         unsigned char ref = line[0];
         if      (ref == 'A' || ref == 'a') ref = 0;
         else if (ref == 'C' || ref == 'c') ref = 1;
         else if (ref == 'G' || ref == 'g') ref = 2;
         else if (ref == 'T' || ref == 't') ref = 3;
         // Convert colors to bases.
         int i = 0;
         for (; i < bytesrd-2; i++) {
            if (line[i+1] == '.') break;
            bases[i] = to_base[(ref = colspace[ref][line[i+1]-48])];
         }
         for (; i < bytesrd-2; i++) {
            bases[i] = 'N';
         }
         // Print converted sequence.
         fprintf(stdout, "%s\n", bases);
      } else if (lineno % div == qline) {
         // Delete first quality score.
         fprintf(stdout, "%s", line+1);
      } else {
         fprintf(stdout, "%s", line);
      }
      lineno++;
   }

   fclose(fin);
   if (call_zcat) {
      close(pipefd[0]);
   }

   return 0;
}
