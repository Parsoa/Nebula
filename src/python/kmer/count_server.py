import io
import os
import pwd
import sys
import copy
import json
import time
import socketserver

from . import (
    bed,
    sets,
    config,
    commons,
    reference,
    counttable,
)

import json
import khmer
import socket
import struct
import colorama
import pybedtools

def colorful_print(*args):
    print(colorama.Fore.CYAN, *args)

# @commons.measure_time
def get_kmer_count(kmer, index):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect(('localhost', 6985 + index))
    s.send(bytearray(kmer, 'ascii'))
    response = s.recv(4) # integer size
    count = struct.unpack('!i', response)[0]
    return count

class CountTableServerHandler(socketserver.BaseRequestHandler):

    sample_counttable = None

    def setup(self):
        return socketserver.BaseRequestHandler.setup(self)

    def handle(self):
        c = config.Configuration()
        buffer = self.request.recv(c.ksize)
        fmt = str(c.ksize) + 's'
        kmer = struct.unpack(fmt, buffer)[0].decode("ascii") 
        count = CountTableServerHandler.sample_counttable.get_kmer_counts(kmer)[0]
        # colorful_print(kmer, ': ', count)
        self.request.send(struct.pack('!i', count))
        return

    def finish(self):
        return socketserver.BaseRequestHandler.finish(self)

class CountTableServer(socketserver.TCPServer):

    def __init__(self, server_address, handler_class):
        socketserver.UnixStreamServer.__init__(self, server_address, handler_class)
        return

    def server_activate(self):
        socketserver.UnixStreamServer.server_activate(self)
        return

    def serve_forever(self):
        while True:
            self.handle_request()
        return

    def handle_request(self):
        return socketserver.TCPServer.handle_request(self)

    def verify_request(self, request, client_address):
        return socketserver.TCPServer.verify_request(self, request, client_address)

    def process_request(self, request, client_address):
        return socketserver.TCPServer.process_request(self, request, client_address)

    def server_close(self):
        return socketserver.TCPServer.server_close(self)

    def finish_request(self, request, client_address):
        return socketserver.TCPServer.finish_request(self, request, client_address)

    def close_request(self, request_address):
        return socketserver.TCPServer.close_request(self, request_address)

if __name__ == '__main__':
    config.configure()
    c = config.Configuration()
    # load k-mer counts
    # this is shared between all children
    colorful_print("loading counttable")
    counttable.export_sample_counttable()
    CountTableServerHandler.sample_counttable = counttable.import_sample_counttable()
    colorful_print("done!")
    #
    children = []
    print('spawning chlidren ...')
    for i in range(0, c.num_threads) :
        pid = os.fork()
        if pid == 0:
            # child
            print('starting server ', i)
            address = ('localhost', 6985 + i)
            server = CountTableServer(server_address = address, handler_class = CountTableServerHandler)
            ip, port = server.server_address
            colorful_print('ip: ', ip, ', port: ', port)
            print('serving forever ', i, ' ...')
            server.serve_forever()
            exit()
        else:
            children.append(pid)
    # only parent process will ever get here
    os.waitpid(children[0], 0)

