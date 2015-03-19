#include "seq_file.h"
SETUP_SEQ_FILE();

int main(int argc, char **argv)
{
  if(argc != 2) exit(EXIT_FAILURE);
  seq_file_t *f = seq_open(argv[1]);
  read_t *r = seq_read_alloc();
  if(f == NULL) exit(EXIT_FAILURE);
  while(seq_read(f,r) > 0)
    printf("%s\t[%lu,%lu,%lu]\n", r->name.b, r->name.end, r->seq.end, r->qual.end);
  seq_close(f);
  seq_read_destroy(r);
  return EXIT_SUCCESS;
}
