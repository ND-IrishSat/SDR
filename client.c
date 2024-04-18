#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <netdb.h>
#include <netinet/in.h>
#include <sys/socket.h>

const char *HOST = "localhost";
const char *PORT = "5050";

int main() {
    /* Lookup server address information */
    struct addrinfo  hints = {
        .ai_family   = AF_UNSPEC,   /* Return IPv4 and IPv6 choices */
        .ai_socktype = SOCK_STREAM, /* Use TCP */
    };
    struct addrinfo *results;
    int status;

    if ((status = getaddrinfo(HOST, PORT, &hints, &results)) != 0) {
        fprintf(stderr, "getaddrinfo failed: %s\n", gai_strerror(status));
        return EXIT_FAILURE;
    }

    /* Create a socket */
    int client_fd = socket(results->ai_family, results->ai_socktype, results->ai_protocol);
    if (client_fd == -1) {
        perror("socket");
        return EXIT_FAILURE;
    }

    /* Connect to the server */
    if (connect(client_fd, results->ai_addr, results->ai_addrlen) == -1) {
        perror("connect");
        close(client_fd);
        return EXIT_FAILURE;
    }
    printf("Hello! Please enter your name:\n");
    char name[BUFSIZ];
    fgets(name, BUFSIZ, stdin);
    send(client_fd, name, strlen(name), 0);
     /* Send floats to the server */
    printf("Enter floats line by line and press Ctrl + D to end: ");
    char input[3*BUFSIZ];
    while(fgets(input, 3*BUFSIZ, stdin) != NULL){
        if (send(client_fd, input, strlen(input), 0) == -1) {
        perror("send");
        close(client_fd);
        return EXIT_FAILURE;
    }
    }

    /* Clean up */
    close(client_fd);
    freeaddrinfo(results);

    return EXIT_SUCCESS;
}
