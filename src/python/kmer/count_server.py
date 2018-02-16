import io
import os
import pwd
import sys
import copy
import json
import time
import argparse
import socketserver

from subprocess import call

from kmer import (
    bed,
    sets,
    config,
    commons,
    counttable,
)

import json
import khmer
import pylru
import socket
import struct
import colorama
import pybedtools

print('importing count_server.py')
# ================================================================================================= #
# Helpers
# ================================================================================================= #

def colorful_print(*args):
    print(colorama.Fore.CYAN, args)

# ================================================================================================= #
# Caching
# ================================================================================================= #

def cache(f):
    cache = pylru.lrucache(config.Configuration.kmer_cache_size)
    hits = 0
    misses = 0
    def wrapper(kmer, index):
        if kmer in cache:
            print('miss')
            nonlocal misses
            misses += 1
            return cache[kmer]
        nonlocal hits
        hits += 1
        print('hit: ', hits / (hits + misses), 'hits: ', hits, ' misses: ', misses)
        cache[kmer] = f(kmer, index) 
        return cache[kmer]
    return wrapper

# ================================================================================================= #
# 
# ================================================================================================= #

def kill():
    call(["pkill", "-f", "count_server"])

def get_kmer_count(kmer, index, ref):
    c = config.Configuration()
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    port = c.count_server_port
    # port = c.reference_count_server_port if ref else c.sample_count_server_port
    s.connect(('localhost', port + index))
    s.send(bytearray(kmer, 'ascii'))
    response = s.recv(4) # integer size
    count = struct.unpack('!i', response)[0]
    return count

# ================================================================================================= #
# 
# ================================================================================================= #

class CountTableServerHandler(socketserver.BaseRequestHandler):

    khmer_counttable = None

    def setup(self):
        return socketserver.BaseRequestHandler.setup(self)

    def handle(self):
        c = config.Configuration()
        buffer = self.request.recv(c.ksize)
        fmt = str(c.ksize) + 's'
        kmer = struct.unpack(fmt, buffer)[0].decode("ascii") 
        count = CountTableServerHandler.khmer_counttable.get_kmer_counts(kmer)[0]
        self.request.send(struct.pack('!i', count))
        return

    def finish(self):
        return socketserver.BaseRequestHandler.finish(self)

class CountTableServer(socketserver.TCPServer):

    # this needs the --counttable arg to work

    @staticmethod
    def run_server():
        c = config.Configuration()
        # 
        counttable.export_counttable()
        CountTableServerHandler.khmer_counttable = counttable.import_counttable()
        # load k-mer counts
        # this is shared between all children
        children = {}
        port = c.count_server_port
        print('spawning chlidren ...')
        for i in range(0, c.max_threads) :
            pid = os.fork()
            if pid == 0:
                # child
                print('starting server ', i)
                address = ('localhost', port + i)
                server = CountTableServer(server_address = address, handler_class = CountTableServerHandler)
                ip, port = server.server_address
                colorful_print('ip: ', ip, ', port: ', port)
                print('serving forever ', i, ' ...')
                server.serve_forever()
                # will never exit
                exit()
            else:
                # parent
                children[pid] = True
        # only parent process will ever get here
        while True:
            (pid, e) = os.wait()
            children.pop(pid, None)
            print(colorama.Fore.RED, 'pid ', pid, 'finished')
            if len(children) == 0:
                break

    def __init__(self, server_address, handler_class):
        socketserver.UnixStreamServer.__init__(self, server_address, handler_class)
        return

    def server_activate(self):
        socketserver.UnixStreamServer.server_activate(self)
        return

    def serve_forever(self):
        while True:
            try:
                self.handle_request()
            except Exception as e:
                print(e)
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

# ================================================================================================= #
#
# ================================================================================================= #

if __name__ == '__main__':
    config.init()
    # 
    CountTableServer.run_server()
