import time
import socket
import struct

from . import (
    config,
)

def colorful_print(*args):
    print(colorama.Fore.CYAN, *args)

if __name__ == '__main__':
    config.configure()

    kmer = 'ATCGATCGATCGATCGATGCCTAGCTAGCTTAGCT'
    while True:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect(('localhost', 6985))
        s.send(bytearray(kmer, 'ascii'))
        response = s.recv(4) # integer size
        count = struct.unpack('!i', response)[0]
        print('count: ', count)
        time.sleep(10)
