'''
A simple interface for network handshaking and socket handling for client/server tasks
'''

import socket
import threading
import errno
import zlib

from socket import error as socket_error

class HandShake():
    ''' 
    A roll-your-own implementation of handshaking loosely based on WebSockets.  HandShake is a parent class 
    containing methods and attributes useful for negotiating common client/server tasks (e.g., file send/recieve, 
    searching, etc...).
    '''
    def __init__(self,*kargs):
        self.HOST = ''
        self.CONNECTION = ''
        self.SEC_KEY = ''
        self.SEC_PROTOCOL = ''
        self.SEC_VERSION = ''
        self.MESSAGE = ''

class ClientConnection(threading.Thread):
    def __init__(self,*kargs):
        self.h = HandShake()
        threading.Thread.__init__(self)
        try:
            self.socket = kargs[0]
            self.connection, self.address = self.socket.accept()
            print " -- client connect from", self.address[0]
        except socket_error as sError:
            raise sError
    def echo(self):
        while True:
            try:
                data = self.connection.recv(256)
            except socket_error as sError:
                raise sError
            if data :
                self.connection.send(data)
    def run(self):
        self.echo()

class GetFile(ClientConnection):
    def __init__(self,*kargs):
        ClientConnection.__init__(self,kargs[0])
    def handshake(self):
        print "hi"

class SendFile(ClientConnection):
    def __init__(self,*kargs):
        ClientConnection.__init__(self,kargs[0])
    def handshake(self):
        print "hi"

class LocateFile(ClientConnection):
    def __init__(self,*kargs):
        ClientConnection.__init__(self,kargs[0])

class Authenticate(ClientConnection):
    def __init__(self,*kargs):
        ClientConnection.__init__(self,kargs[0])

class Server:
    def __init__(self,*kargs):
        self.socket = socket.socket()
        if(len(kargs)>0):
            self.port = kargs[0]
            self.maxConnect = kargs[1] # sets the number of queued clients on our single-thread stack while we work with our active connection
        else:
            self.port = 65000
            self.maxConnect = 5
        try:
            self.socket.bind(('',self.port))  # bind to localhost on port N
        except socket_error as sError:
            raise sError
    def listen(self):
        try:
            self.socket.listen(self.maxConnect)
        except socket_error as sError:
            raise sError
        self.client_list = []
        while True:
            self.client_list.append(ClientConnection(self.socket).start())
    def close(self):
        self.socket.close()

if __name__ == "__main__":
    s=Server(12345,10)
    s.listen()
