#include "seq_file.h"
SETUP_SEQ_FILE();

#define print_prompt() printf(">"); fflush(stdout)

int main(int argc, char **argv)
{
  (void)argv;

  if(argc != 2) {
    fprintf(stderr, "usage: seq_cmdline\n");
    exit(EXIT_FAILURE);
  }

  seq_file_t *file = seq_open_fh(stdin,0,0);

  if(file == NULL)
  {
    fprintf(stderr, "Error: couldn't open stdin\n");
    exit(EXIT_FAILURE);
  }

  read_t *read = seq_read_alloc();

  print_prompt();
  while(seq_read(file, read) > 0)
  {
    printf("%s\n", read->seq.b);
    print_prompt();
  }

  seq_close(file);
  seq_read_destroy(read);

  return EXIT_SUCCESS;
}
