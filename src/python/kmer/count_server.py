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

def get_kmer_count(kmer):
    print('kire khar')

class CountTableServerHandler(socketserver.BaseRequestHandler):

    def setup(self):
        colorful_print("loading counttable")
        counttable.export_sample_counttable()
        self.counttable = counttable.import_sample_counttable()
        colorful_print("done!")
        return socketserver.BaseRequestHandler.setup(self)

    def handle(self):
        colorful_print("handling request")
        data = self.request.recv(1024)
        colorful_print('recv: ', data)
        # count = self.counttable.get_kmer_counts(data)
        # colorful_print('count: ', count)
        self.request.send(data)
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
    address = ('localhost', 6985)
    server = CountTableServer(server_address = address, handler_class = CountTableServerHandler)
    ip, port = server.server_address
    colorful_print('ip: ', ip, ', port: ', port)
    server.serve_forever()
