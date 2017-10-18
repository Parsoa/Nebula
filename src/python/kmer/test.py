import time
import socket

if __name__ == '__main__':

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect(('localhost', 6985))
    while True:
        time.sleep(2)
        len_sent = s.send(b'hello kire khar')
        response = s.recv(len_sent)