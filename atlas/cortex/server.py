#!/usr/bin/env python

## Adapted from https://github.com/mcveanlab/mccortex/blob/master/scripts/mccortex-server.py

from __future__ import print_function

import os
import sys
import time
import signal
import requests
from Bio.Seq import Seq

from subprocess import Popen, PIPE

# These should work on python2 after 'pip install --user future'
from http.server import BaseHTTPRequestHandler
import socketserver
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

#
# Start server on port 2306, with link and graph files:
#   python mccortex-server.py 2306 --coverages --edges -p l.ctp.gz a.ctx b.ctx
# Query a kmer:
#   curl localhost:2096/CAGTGGCCA
# Response:
#   { "key": "CAGTGGCCA", "colours": [1], "left": "T", "right": "T", "edges": "88", "links": [] }
#

def check_mccortex_alive(proc):
    if proc.poll() is not None:
        logger.info("McCortex quit [%i]" % proc.returncode, file=sys.stdout)
        sys.exit(1)

def query_mccortex(proc,kmer):
    print(kmer, file=proc.stdin)
    proc.stdin.flush()
    check_mccortex_alive(proc)
    line = proc.stdout.readline()
    check_mccortex_alive(proc)
    return line

# when we start mccortex we set it to ignore interrupt signal, we handle it.
def preexec_function():
    # Ignore the SIGINT signal by setting the handler to the standard
    # signal handler SIG_IGN.
    signal.signal(signal.SIGINT, signal.SIG_IGN)

class WebServer(object):

    def __init__(self, port, args):
        self.port = port
        self.args = args
        self.mccortex = None
        

    def start(self):
        self.mccortex = self._start_mccortex()

    def serve(self):
        mccortex = self.mccortex
        class McCortexHTTPServer(BaseHTTPRequestHandler):

            def do_GET(self):
                self.send_response(200)
                self.send_header("Content-type", "text/html")
                self.end_headers()
                if len(self.path) < 4 or len(self.path) > 300:
                    jsonstr = "{\"error\": \"Webserver: bad query\"}\n"
                else:
                    jsonstr = query_mccortex(mccortex, self.path[1:])
                self.wfile.write(jsonstr.encode("UTF-8"))        
        try:
            self.httpd = socketserver.TCPServer( ("", self.port), McCortexHTTPServer)
        except Exception as e:
            logger.error("Cannot start HTTP server: %s" % str(e))
            sys.exit()           
        try:
            self.httpd.serve_forever()
        except:
            logger.info("Exited server") 

    def stop(self):
        self.httpd.server_close()
        self._stop_mccortex()
        logger.info("Closed succesfully.")   

    def _start_mccortex(self):
        script_dir = os.path.dirname(os.path.realpath(__file__))
        # Adding two lists together appends one to the other
        try:
            proc = Popen(["mccortex31", "server", "--single-line", "-q"] + self.args,
                          stdin=PIPE, stdout=PIPE, universal_newlines=True,
                          preexec_fn = preexec_function)
        except Exception as e:
            logger.error("Couldn't start McCortex: %s " % str(e))
            sys.exit(1)

        # Give test query to check it works
        check_mccortex_alive(proc)
        resp = query_mccortex(proc, "hi")
        return proc

    def _stop_mccortex(self):
        print("q\n", file=self.mccortex.stdin)
        self.mccortex.stdin.close()
        # sleep until process has closed
        while self.mccortex.poll() is None:
            time.sleep(1)
        logger.info("McCortex exited with: %s " % str(self.mccortex.poll()))            

class McCortexQuery(object):

    def __init__(self, port, seeds = []):
          self.base_url = "http://localhost:%i/" % port
          self.seeds = seeds

    def query(self, kmer):
        return McCortexQueryResult(kmer, requests.get(self.base_url + kmer).json(), seeds = self.seeds)

class McCortexQueryResult(object):

    def __init__(self, query_kmer, data, seeds = []):
      self.query_kmer = query_kmer
      self.data = data
      self.seeds = seeds

    def __str__(self):
      return str(self.data)

    def __repr__(self):
      return str(self.data)

    def forward(self):
      forward_kmers = []
      if self.complement:
          for l in self.left:
            forward_kmers.append(str(Seq(l + self.data["key"][:-1]).reverse_complement()))
      else:
          for r in self.right:
            forward_kmers.append(str(self.data["key"][1:] + r))
      if len(forward_kmers) > 1:
          forward_kmers = [f for f in forward_kmers if f in self.seeds]
      return forward_kmers

    @property
    def right(self):
      return self.data["right"]

    @property
    def left(self):
      return self.data["left"]     

    @property
    def kmer(self):
      return self.data["key"]

    @property
    def depth(self):
      try:
          return int(self.data["edges"])
      except ValueError:
          try:
              return int(self.data["edges"][-1])
          except ValueError:
              return 0
    @property
    def complement(self):
      return self.kmer != self.query_kmer

# def main():
#     if len(sys.argv) < 3 or not sys.argv[1].isdigit():
#         print("usage: %s <port> [mccortex args]" % (sys.argv[0]))
#         sys.exit(-1)

#     port = int(sys.argv[1])
#     args = sys.argv[2:]

#     start_web_server(port, args)

# if __name__ == '__main__':
#     main()