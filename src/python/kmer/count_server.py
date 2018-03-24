import io
import os
import pwd
import sys
import copy
import json
import time
import argparse
import SocketServer as socketserver

from subprocess import call

from kmer import (
    config,
    counttable,
)

#import pylru
import socket
import struct

print('importing count_server.py')
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

    counts_provider = None

    def setup(self):
        return socketserver.BaseRequestHandler.setup(self)

    def handle(self):
        c = config.Configuration()
        buffer = self.request.recv(c.ksize)
        fmt = str(c.ksize) + 's'
        kmer = struct.unpack(fmt, buffer)[0].decode("ascii") 
        count = CountTableServerHandler.counts_provider.get_kmer_count(kmer)
        self.request.send(struct.pack('!i', count))
        return

    def finish(self):
        return socketserver.BaseRequestHandler.finish(self)

class CountTableServer(socketserver.TCPServer):

    def __init__(self, server_address, handler_class):
        socketserver.TCPServer.__init__(self, server_address, handler_class)

    def server_activate(self):
        socketserver.TCPServer.server_activate(self)
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

def run_server():
    c = config.Configuration()
    # 
    if c.is_dummy:
        CountTableServerHandler.counts_provider = counttable.DummyCountsProvider()
    #elif c.khmer:
    #    CountTableServerHandler.counts_provider = counttable.KhmerCountsProvider()
    elif c.jellyfish:
        CountTableServerHandler.counts_provider = counttable.JellyfishCountsProvider()
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
            server = CountTableServer(address, CountTableServerHandler)
            ip, port = server.server_address
            print('ip: ', ip, ', port: ', port)
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

# ================================================================================================= #
#
# ================================================================================================= #

if __name__ == '__main__':
    config.init()
    run_server()
