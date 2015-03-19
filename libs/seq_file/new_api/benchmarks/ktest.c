#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
  if (argc == 1) return 1;
  gzFile fp = gzopen(argv[1], "r");
  kseq_t *ks = kseq_init(fp);
  while (kseq_read(ks) >= 0)
    printf("%s\t[%lu,%lu,%lu]\n", ks->name.s, ks->name.l, ks->seq.l, ks->qual.l);
  kseq_destroy(ks);
  gzclose(fp);
  return 0;
}
