#include "buffered_input.h"

BUFFERED_INPUT_SETUP();

int main(int argc, char **argv)
{
  if(argc != 2) exit(EXIT_FAILURE);

  buffer_t *in = buffer_alloc(10);
  buffer_t *buf = buffer_alloc(10);
  // gzFile f = gzopen(argv[1],"r");
  FILE *f = fopen(argv[1],"r");
  if(f == NULL) exit(EXIT_FAILURE);

  int i;
  for(i = 0; freadline_buf(f,in,buf) > 0; i++)
  {
    printf("line %i: '%s'\n", i, buf->b);
  }

  fclose(f);
  buffer_destroy(buf);
}
