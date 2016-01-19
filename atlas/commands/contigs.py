from __future__ import print_function
## Using seed kmers from a gene panel return a DFS through the graph
import sys
sys.path.append('/home/phelimb/git/atlas-core')
from atlas.cortex.server import WebServer
from atlas.cortex.server import McCortexQuery
import socket
import threading
import logging
logger = logging.getLogger(__name__)

# def run(parser, args):
#     args = parser.parse_args()
#     check_args(args)  
#     Genotyper(args).run()

def get_open_port():
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(("",0))
        s.listen(1)
        port = s.getsockname()[1]
        s.close()
        return port




port =  (get_open_port())
logger.info("Running on port %i " % port)
wb = WebServer(port, args = [ "/data2/users/phelim/data/gn/atlas/bins/k31/C00026396.ctx"] )
logger.info("Loading binary")
wb.start()
## Serve on a thread
logger.info("Starting sever")
server = threading.Thread(target=wb.serve)
server.start()

logger.info("Walking the graph")
query = McCortexQuery(port = port)
print (query.query("CTCCACCACTTGCCCCTGCTTCATCACCAGC"))
logger.info("Cleaning up")
wb.stop()




