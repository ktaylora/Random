"""
A simple interface for network handshaking and socket handling for client/server tasks
"""

import socket
import threading
import errno
import zlib

from socket import error as socket_error


class Handshake():
    def __init__(self, *args):
        """
        A roll-your-own implementation of handshaking loosely based on WebSockets.  HandShake is a parent class
        containing methods and attributes useful for negotiating common client/server tasks (e.g., file send/recv,
        searching, etc...).
        :rtype: object
        """
        self.HOST = ''
        self.CONNECTION = ''
        self.SEC_KEY = ''
        self.SEC_PROTOCOL = ''
        self.SEC_VERSION = ''
        self.MESSAGE = ''

    def processHandshake(self):
        known = False
        if self.MESSAGE == "001":
            known = True
            print "hi"
        elif self.MESSAGE == "002":
            known = True
            print "hi"


class ClientConnection(threading.Thread):
    def __init__(self, *args):
        """
        parent class providing exposure to handshaking and threading for client-handling tasks
        :param args:
        :return:
        """
        threading.Thread.__init__(self)
        try:
            self.socket = args[0]
            self.connection, self.address = self.socket.accept()
            print " -- client connect from", self.address[0]
        except socket_error as sError:
            raise sError

    def run(self):
        """
        default run method provides a simple threaded echo server.  Default method should be overwritten by derived
        children to provide more complex functionality
        :return NULL:
        """
        while True:
            try:
                data = self.connection.recv(256)
            except socket_error as sError:
                raise sError
            if data:
                self.connection.send(data)


class ServerConnection(threading.Thread):
    def __init__(self, *args):
        """
        parent class providing exposure to handshaking and threading for client-handling tasks
        :param args:
        :return:
        """
        self.h = Handshake()
        threading.Thread.__init__(self)


class SGetFile(ClientConnection):
    def __init__(self, *args):
        """
        GetFile parent class serves/requests files to a connected client upon request
        :param args:
        :return:
        """
        ClientConnection.__init__(self, args[0])


class CGetFile(ServerConnection):
    def __init__(self, *args):
        """
        GetFile parent class serves/requests files to a connected client upon request
        :param args:
        :return:
        """
        ServerConnection.__init__(self, args[0])


class SLocateFile(ClientConnection):
    def __init__(self, *args):
        """
        client-related LocateFile interface that accepts a filename and attempts to locate a file on the server-side,
        using UNIX tools.
        :param args:
        :return:
        """
        ClientConnection.__init__(self, args[0])


class SAuthenticate(ClientConnection):
    def __init__(self, *args):
        """
        client-related Authenicate interface that handles a user request for authentication
        :param args:
        :return:
        """
        ClientConnection.__init__(self, args[0])


class CAuthenticate(ServerConnection):
    def __init__(self, *args):
        """
        client-related Authenicate interface that handles a user request for authentication
        :param args:
        :return:
        """
        ServerConnection.__init__(self, args[0])


class SGetClientHandshake(ClientConnection, Handshake):
    def __init__(self, *args):
        Handshake.__init__(self)
        ClientConnection.__init__(self, args[0])

    def run(self):
        """
        default run method provides a simple threaded echo server.  Default method should be overwritten by derived
        child classes to provide more complex functionality
        :return NULL:
        """
        while True:
            try:
                data = self.connection.recv(256)
            except socket_error as sError:
                raise sError
            if data:
                self.connection.send(data)


class Server:
    def __init__(self, *args):
        """

        :param args:
        :return:
        """
        self.socket = socket.socket()
        if len(args) > 0:
            self.port = args[0]
            self.maxConnect = args[
                1]  # sets the number of queued clients on our single-thread stack while we work with our active connection
        else:
            self.port = 65000
            self.maxConnect = 5
        try:
            self.socket.bind(('', self.port))  # bind to localhost on port N
        except socket_error as sError:
            raise sError

    def listen(self):
        """

        :rtype: object
        """
        try:
            self.socket.listen(self.maxConnect)
        except socket_error as sError:
            raise sError

        client_list = []
        while True:
            client_list.append(SGetClientHandshake(self.socket).start())

    def close(self):
        self.socket.close()


if __name__ == "__main__":
    try:
        s = Server(12345, 10)
        s.listen()
    except KeyboardInterrupt as e:
        s.close()
    except Exception as e:
        print e
