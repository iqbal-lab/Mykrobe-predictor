
#include "seq_file.h"

SETUP_SEQ_FILE();

int main(int argc, char **argv)
{
  if(argc != 2)
  {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  int min, max;
  char get_quals = seq_get_qual_limits(argv[1], 200, &min, &max);
  if(get_quals > 0) printf("qual range %i - %i\n", min, max);
  if(get_quals == 0) printf("no qual values in file\n");
  if(get_quals < 0) printf("IO Error: cannot read qual values\n");

  read_t *r = seq_read_alloc();
  seq_file_t *f = seq_open(argv[1]);

  if(f == NULL)
  {
    fprintf(stderr, "Cannot open file\n");
    exit(EXIT_FAILURE);
  }

  int s = seq_read(f,r);
  if(s > 0)
  {
    if(seq_is_sam(f)) printf("file is sam\n");
    if(seq_is_bam(f)) printf("file is bam\n");
    if(seq_is_fasta(f)) printf("file is fasta\n");
    if(seq_is_fastq(f)) printf("file is fastq\n");
    if(seq_is_plain(f)) printf("file is plain\n");
    if(seq_use_gzip(f)) printf("file is using gzip\n");
  }
  else if (s < 0) printf("Error occurred reading file\n");
  else printf("Cannot get any reads from file: %s\n", argv[0]);

  seq_close(f);
  seq_read_destroy(r);

  return EXIT_SUCCESS;
}
