from __future__ import print_function
import os
import sys
import json
import time
import signal
from Bio.Seq import Seq
from atlas.utils import median
from subprocess import Popen, PIPE
from http.server import BaseHTTPRequestHandler
import socketserver
import logging
from pprint import pprint
import copy
logger = logging.getLogger(__name__)


# Adapted from
# https://github.com/mcveanlab/mccortex/blob/master/scripts/mccortex-server.py
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


def query_mccortex(proc, kmer):
    print(kmer, file=proc.stdin)
    proc.stdin.flush()
    check_mccortex_alive(proc)
    line = proc.stdout.readline()
    # Trim off prompt text
    if line[0:2] == "> ":
        line = line[2:len(line)]
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
        self.httpd = None

    def start(self):
        self.mccortex = self._start_mccortex()

    def serve(self):
        mccortex = self.mccortex
        logger.debug("Starting server")

        class McCortexHTTPServer(BaseHTTPRequestHandler):

            def log_message(self, format, *args):
                # Remove logging
                pass

            def do_GET(self):
                self.send_response(200)
                self.send_header("Content-type", "text/html")
                self.end_headers()
                if len(self.path) < 4 or len(self.path) > 300:
                    jsonstr = "{\"error\": \"Webserver: bad query\"}\n"
                else:
                    jsonstr = query_mccortex(mccortex, self.path[1:])
                try:
                    self.wfile.write(jsonstr.encode("UTF-8"))
                except UnicodeDecodeError:
                    self.wfile.write("{}".encode("UTF-8"))
        try:
            self.httpd = socketserver.TCPServer(
                ("", self.port), McCortexHTTPServer)
        except Exception as e:
            logger.error("Cannot start HTTP server: %s" % str(e))
            sys.exit()
        try:
            self.httpd.serve_forever()
        except:
            logger.info("Exited server")

    def stop(self):
        if self.httpd:
            self.httpd.server_close()
        self._stop_mccortex()
        logger.info("Closed succesfully.")

    def _start_mccortex(self):
        script_dir = os.path.dirname(os.path.realpath(__file__))
        # Adding two lists together appends one to the other
        try:
            proc = Popen(["mccortex31",
                          "server",
                          "--single-line",
                          "--coverages"] + self.args,
                         stdin=PIPE,
                         stdout=PIPE,
                         universal_newlines=True,
                         preexec_fn=preexec_function)
        except Exception as e:
            logger.error(
                "Couldn't start McCortex is mccortex31 in path? : %s " %
                str(e))
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

    def __init__(self, proc):
        # self.base_url = "http://localhost:%i/" % port
        self.proc = proc

    def query(self, kmer, known_kmers=[]):
        return McCortexQueryResult(
            kmer,
            json.loads(
                query_mccortex(
                    self.proc,
                    kmer)),
            known_kmers=known_kmers)


class McCortexQueryResult(object):

    def __init__(self, query_kmer, data, known_kmers=[]):
        self.query_kmer = query_kmer
        self.data = data
        self.known_kmers = known_kmers

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return str(self.data)

    def forward(self, suggested_kmer=None):
        forward_kmers = []
        if self.complement:
            for l in self.left:
                forward_kmers.append(
                    str(Seq(l + self.data["key"][:-1]).reverse_complement()))
        else:
            for r in self.right:
                forward_kmers.append(str(self.data["key"][1:] + r))
        if suggested_kmer is not None and suggested_kmer in forward_kmers:
            return [suggested_kmer]
        if len(forward_kmers) > 1 and self.known_kmers:
            _forward_kmers = [
                f for f in forward_kmers if f in self.known_kmers]
            if len(_forward_kmers) > 0:
                forward_kmers = _forward_kmers
            assert len(forward_kmers) > 0
        return forward_kmers

    def reverse(self, suggested_kmer=None):
        reverse_kmers = []
        if self.complement:
            for l in self.right:
                reverse_kmers.append(
                    str(Seq(self.data["key"][1:] + l).reverse_complement()))
        else:
            for r in self.left:
                reverse_kmers.append(r + str(self.data["key"][:-1]))
        if suggested_kmer is not None and suggested_kmer in reverse_kmers:
            return [suggested_kmer]
        if len(reverse_kmers) > 1 and self.known_kmers:
            _reverse_kmers = [
                f for f in reverse_kmers if f in self.known_kmers]
            if len(_reverse_kmers) > 0:
                reverse_kmers = _reverse_kmers
            assert len(reverse_kmers) > 0
        return reverse_kmers

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
            return int(self.data["colours"][0])
        except ValueError:
            logger.error(
                "parsing edges %s:%s " %
                (self.data["key"], self.data["edges"]))
            return None

    @property
    def complement(self):
        return self.kmer != self.query_kmer


class GraphWalker(object):

    def __init__(self, proc, kmer_size=31, print_depths=False):
        self.mcq = McCortexQuery(proc)
        self.queries = {}
        self.kmer_size = kmer_size
        self.print_depths = print_depths

    def _count_k(self, k, count):
        try:
            count[k] += 1
        except KeyError:
            count[k] = 1
        return count

    def _path_ended(self, paths, i):
        return paths[i]["dna"][-1] == "*"

    def _get_next_kmer_right(self, paths, i):
        return paths[i]["dna"][-self.kmer_size:]

    def _get_next_kmer_left(self, paths, i):
        return paths[i]["dna"][:self.kmer_size]

    def _reached_end_of_path(self, k, end_kmers):
        return k in end_kmers

    def _mark_path_and_ended(self, paths, i):
        paths[i]["dna"] = paths[i]["dna"] + "*"
        return paths

    def _make_query(self, k, known_kmers=[]):
        try:
            q = self.queries[k]
        except KeyError:
            try:
                q = self.mcq.query(k, known_kmers=known_kmers)
                self.queries[k] = q
            except ValueError as e:
                q = None
                logging.error(str(e))
                self.queries[k] = None
        return q

    def _query_is_valid(self, q):
        return q is not None and q.data.get("key")

    def _get_next_kmers_right(
            self,
            i,
            q,
            k,
            paths,
            repeat_kmers,
            known_kmers,
            count):
        kmers = q.forward(suggested_kmer=repeat_kmers.get(k, {}).get(count[k]))
        if q.depth is not None:
            if self.print_depths:
                paths[i]["depth"].append(q.depth)
            if len(kmers) > 1:
                depth = max(q.depth * 0.1, 10)
                kmers = [
                    k for k in kmers if self.mcq.query(
                        k,
                        known_kmers=known_kmers).depth > depth]
        if len(kmers) > 1:
            kmers = [k for k in kmers if k not in paths[i]["dna"]]
        return kmers

    def _get_next_kmers_left(
            self,
            i,
            q,
            k,
            paths,
            repeat_kmers,
            known_kmers=[],
            count={}):
        kmers = q.reverse()
        if q.depth is not None:
            if len(kmers) > 1:
                depth = max(q.depth * 0.1, 10)
                kmers = [
                    k for k in kmers if self.mcq.query(
                        k,
                        known_kmers=known_kmers).depth > depth]
        if len(kmers) > 1:
            kmers = [k for k in kmers if k not in paths[i]["dna"]]
        return kmers

    def _check_if_next_kmers_are_valid(self, kmers, paths, i):
        if len(kmers) < 1:
            paths = self._mark_path_and_ended(paths, i)
        return paths

    def _split_paths(self, i, kmers, paths):
        for j, kmer in enumerate(kmers):
            if j == 0:
                paths[i]["dna"] = paths[i]["dna"] + kmer[-1]
            else:
                paths[sorted(
                    paths.keys())[-j]]["dna"] = paths[sorted(paths.keys())[-j]]["dna"] + kmer[-1]
        return paths

    def _init_paths(self, seed, N_left):
        if self.print_depths:
            return {0: {"N_left": N_left,
                        "dna": seed[:self.kmer_size],
                        "start_kmer": seed[:self.kmer_size],
                        "depth": []}}
        else:
            return {0: {"N_left": N_left,
                        "dna": seed[:self.kmer_size],
                        "start_kmer": seed[:self.kmer_size]}}

    def breath_first_search(self, N, seed, end_kmers=[],
                            known_kmers=[], repeat_kmers={},
                            N_left=0):
        paths = self._init_paths(seed, N_left)
        count = {}
        for _ in range(N_left):
            k = self._get_next_kmer_left(paths, 0)
            q = self._make_query(k)
            if self._query_is_valid(q):
                kmers = self._get_next_kmers_left(
                    0,
                    q,
                    k,
                    paths,
                    repeat_kmers,
                    known_kmers,
                    count)
                if not len(kmers) == 1:
                    break
                else:
                    # Add left
                    paths[0]["dna"] = kmers[0][0] + paths[0]["dna"]
        #
        #
        for _ in range(N - N_left):
            for i in paths.keys():
                if not self._path_ended(paths, i):
                    k = self._get_next_kmer_right(paths, i)
                    count = self._count_k(k, count)
                    if self._reached_end_of_path(k, end_kmers):
                        paths = self._mark_path_and_ended(paths, i)
                    else:
                        q = self._make_query(k, known_kmers)
                        if self._query_is_valid(q):
                            kmers = self._get_next_kmers_right(
                                i,
                                q,
                                k,
                                paths,
                                repeat_kmers,
                                known_kmers,
                                count)
                            paths = self._check_if_next_kmers_are_valid(
                                kmers,
                                paths,
                                i)
                            # Create new paths
                            paths = self.create_new_paths(
                                paths,
                                i,
                                kmers,
                                origin=q.data.get("key"))
                            paths = self._split_paths(i, kmers, paths)
                        else:
                            paths = self._mark_path_and_ended(paths, i)

            # DEBUG prints progress every N bps
            # if i % 250 == 0:
            #     with open("/tmp/paths.json", "w") as outfile:
            #         json.dump(paths, outfile)
            ###

        keep_paths = {}
        for k, v in paths.items():
            v["len_dna"] = len(v["dna"])
            v["prot"] = str(Seq(v["dna"].rstrip("*")).translate(11))
            v["len_prot"] = len(v["prot"])
            if self.print_depths:
                paths[i]["median_depth"] = median(paths[i]["depth"])
                paths[i]["min_non_zero_depth"] = min(paths[i]["depth"])
                paths[i]["depth"] = "-".join([str(x)
                                              for x in paths[i]["depth"]])

            if v["len_dna"] == N + 1 and v["dna"][-1] == "*":
                keep_paths[k] = v
        return keep_paths.values()

    def create_new_paths(self, paths, i, kmers, origin):
        num = len(kmers) - 1
        num_paths = max(paths.keys())
        if num > 0:
            logger.debug("Branch point")
            logger.debug("Origin %s" % origin)
            logger.debug("Options %s" % ",".join(kmers))
        for j in range(num):
            paths[num_paths + j + 1] = copy.copy(paths[i])
            paths[num_paths + j + 1]["start_kmer"] = kmers[1]
            if kmers[1] in [d["start_kmer"] for d in paths.values()]:
                print (kmers, paths)
                raise ValueError("Going around in circles?")
        return paths


# def main():
#     if len(sys.argv) < 3 or not sys.argv[1].isdigit():
#         print("usage: %s <port> [mccortex args]" % (sys.argv[0]))
#         sys.exit(-1)

#     port = int(sys.argv[1])
#     args = sys.argv[2:]

#     start_web_server(port, args)

# if __name__ == '__main__':
#     main()
