/* server.c: simple TCP echo server  */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netdb.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

const char *HOST = NULL;
const char *PORT = "5050";

int main(int argc, char *argv[]) {
    /* Lookup server address information */
    struct addrinfo  hints = {
        .ai_family   = AF_UNSPEC,   /* Return IPv4 and IPv6 choices */
        .ai_socktype = SOCK_STREAM, /* Use TCP */
        .ai_flags    = AI_PASSIVE,  /* Use all interfaces */
    };
    struct addrinfo *results;
    int status;
	/* use getaddrinfo to lookup server address info based on host and port*/
    if ((status = getaddrinfo(HOST, PORT, &hints, &results)) != 0) {
    	fprintf(stderr, "getaddrinfo failed: %s\n", gai_strerror(status));
	return EXIT_FAILURE;
    }
    /* For each server entry, allocate socket and try to connect */
    int server_fd = -1;
    for (struct addrinfo *p = results; p != NULL && server_fd < 0; p = p->ai_next) {
	/* Allocate socket */
	if ((server_fd = socket(p->ai_family, p->ai_socktype, p->ai_protocol)) < 0) {
	    fprintf(stderr, "Unable to make socket: %s\n", strerror(errno));
	    continue;
	}

	/* Bind socket */
	if (bind(server_fd, p->ai_addr, p->ai_addrlen) < 0) {
	    fprintf(stderr, "Unable to bind: %s\n", strerror(errno));
	    close(server_fd);
	    server_fd = -1;
	    continue;
	}
	
    	/* Listen to socket */
	if (listen(server_fd, SOMAXCONN) < 0) {
	    fprintf(stderr, "Unable to listen: %s\n", strerror(errno));
	    close(server_fd);
	    server_fd = -1;
	    continue;
	}
    }

    /* Release allocate address information */
    freeaddrinfo(results);

    if (server_fd < 0) {
		printf("failure\n");
    	return EXIT_FAILURE;
    }

	//puts("about to accept");
    /* Process incoming connections */
	/* accepting incoming connectuon using the accept() function, when a client connects it opens a file stream from the socket file descriptor */
    while (1) {
		struct sockaddr client_addr;
    	socklen_t client_len = sizeof(struct sockaddr);

        /* Accept incoming connection */
    	int client_fd = accept(server_fd, &client_addr, &client_len);
    	if (client_fd < 0) {
    	    fprintf(stderr, "Unable to accept: %s\n", strerror(errno));
    	    continue;

		break;
	}

	puts("have a client");

    /* Open file stream from socket file descriptor */
	FILE *client_file = fdopen(client_fd, "w+");
	if (!client_file) {
    	    fprintf(stderr, "Unable to fdopen: %s\n", strerror(errno));
    	    close(client_fd);
    	    continue;
	}
	//printf("after open file stream");
		/* Read the number of bytes sent by the client */
		char name[BUFSIZ];
		fgets(name, BUFSIZ, client_file);
		name[strlen(name) - 1] ='\0';
		printf("Hello %s!\n", name);
		//fputs(name, stdout);
		char buffer[3*BUFSIZ];
		/*  constantly read a number */
		while(fgets(buffer, 3*BUFSIZ, client_file)){
			fputs(buffer, stdout);
		}
        /* Close connection */
        fclose(client_file);
		close(client_fd);
    }

    return EXIT_SUCCESS;
}

/* vim: set expandtab sts=4 sw=4 ts=8 ft=c: */
