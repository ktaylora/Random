'''
A simple interface for network handshaking and socket handling for client/server tasks
'''

import socket
import getpass

class Keyfile:
    def __init__(self):

class Client:
    def __init__(self):

class Server:
    def __init__(self,*kargs):
        self.socket = socket.socket() # by default, we will assume IP4 SOCK_STREAM
        self.clientConn = ""
        self.address = ""
        if(len(kargs)>0):
            self.port = kargs[0]
            self.maxConnect = kargs[1] # sets the number of queued clients on our single-thread stack while we work with our active connection
        else:
            self.port = 65000
            self.maxConnect = 5

        try:
            self.socket.bind(('',self.port))  # bind to localhost on port N
        except Exception as e:
            print e

    def listen(self):
        self.socket.listen(self.maxConnect)
        while True:
            self.clientConn, self.address = self.socket.accept()
            print " -- client connect from", self.address[0]
            data = self.clientConn.recv(200)
            if data:
                print " -- input recv. from", self.address[0]
                self.clientConn.send(data)
            self.clientConn.close()

    def close(self):
        self.socket.close()
