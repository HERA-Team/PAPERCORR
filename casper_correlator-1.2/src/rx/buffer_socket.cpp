#include <signal.h>
#include <syslog.h>
#include "include/buffer_socket.h"

int default_callback(char *data, size_t size, void *userdata) {
    printf("    Readout packet of size %lu\n", size);
    return 0;
}

void init_buffer_socket(BufferSocket *bs, size_t item_count, size_t payload_len) {
    // Initialize a BufferSocket
    bs->buf = ring_buffer_create(item_count, payload_len);
    set_callback(bs, &default_callback);
    bs->run_threads = 0;
    bs->userdata = NULL;
}

void free_buffer_socket(BufferSocket *bs) {
    // Free all memory allocated for a BufferSocket
    stop(bs);
    ring_buffer_delete(bs->buf);
}

void set_callback(BufferSocket *bs, int (*cb_func)(char *, size_t, void *)) {
    /* Set a callback function for handling data out of ring buffer */
    bs->callback = cb_func;
}

int start(BufferSocket *bs, int port) {
    /* Start socket => buffer and buffer => callback threads */
    if (bs->run_threads != 0) {
        fprintf(stderr, "BufferSocket already running.");
        return 1;
    }
    bs->port = port;
    bs->run_threads = 1;
    pthread_create(&bs->net_thread, NULL, net_thread_function, bs);
    pthread_create(&bs->data_thread, NULL, data_thread_function, bs);
    return 0;
}

int stop(BufferSocket *bs) {
    /* Send halt signal for net/data threads, then join them */
    if (!bs->run_threads) return 1;
    bs->run_threads = 0;
    pthread_join(bs->net_thread, NULL);
    pthread_join(bs->data_thread, NULL);
    return 0;
}
    

void *data_thread_function(void *arg) {
    /* This thread reads data out of a ring buffer through a callback */
    BufferSocket *bs = (BufferSocket *)arg;
    RING_ITEM *this_slot;
    struct timespec ts;

    while (bs->run_threads) {
        this_slot = bs->buf->read_ptr;
        if (clock_gettime(CLOCK_REALTIME, &ts) == -1) {
            fprintf(stderr, "Data: clock_gettime returned nonzero.\n");
            bs->run_threads = 0;
            continue;
        }
        ts.tv_nsec += 10000000;     // 10 ms
        // Wait for next buffer slot to fill up
        if (sem_timedwait(&this_slot->read_mutex, &ts) == -1) continue;
         //printf("Reading in a packet: size=%d slot=%d\n", this_slot->size, this_slot - bs->buf->list_ptr);
        // Feed data from buffer slot to callback function
        if (bs->callback((char *)this_slot->data, 
                this_slot->size, bs->userdata) != 0) {
            fprintf(stderr, "Data: Callback returned nonzero.\n");
            //bs->run_threads = 0;
        } else {
            // Release this slot for writing
            sem_post(&this_slot->write_mutex);
            bs->buf->read_ptr = this_slot->next;
        }
    }
    return NULL;
}

// TIMEOUT_USEC == 10*1000 us == 10 ms
#define TIMEOUT_USEC (10*1000)

void *net_thread_function(void *arg) {
    /* This thread puts data into a ring buffer from a socket*/
    BufferSocket *bs = (BufferSocket *)arg;
    RING_ITEM *this_slot;

    socket_t sock = setup_network_listener((short) bs->port);
    SA_in addr; // packet source's address
    socklen_t addr_len = sizeof(addr);
    ssize_t num_bytes = 0;
    fd_set readset;
    struct timeval tv;
    struct timespec ts;
    long timeouts = 0;

    // If sock open fails, end all threads
    if (sock == -1) {
        fprintf(stderr, "Unable to open socket.\n");
        bs->run_threads = 0;
        return NULL;
    }

    while (bs->run_threads) {
        this_slot = bs->buf->write_ptr;
        if (clock_gettime(CLOCK_REALTIME, &ts) == -1) {
            fprintf(stderr, "Net: clock_gettime returned nonzero.\n");
            bs->run_threads = 0;
            continue;
        }
        ts.tv_nsec += TIMEOUT_USEC*1000;
        // Wait for next buffer slot to open up
        if (sem_timedwait(&this_slot->write_mutex, &ts) == -1) continue;
        
        // Poll to see if socket has data
        while (bs->run_threads) {
            FD_ZERO(&readset);
            FD_SET(sock, &readset);
            tv.tv_sec = 0; tv.tv_usec = TIMEOUT_USEC;
            num_bytes = select(sock + 1, &readset, NULL, NULL, &tv);
            // Read data from socket into ring buffer
            if (num_bytes > 0) {
                num_bytes = recvfrom(sock, this_slot->data, 
                    bs->buf->buffer_size / bs->buf->list_length, 
                    0, (SA *)&addr, &addr_len);
                timeouts = 0;
                break;
            } else if (num_bytes < 0) {
                if (errno == EINTR) continue;
                fprintf(stderr, "Unable to receive packets.\n");
                bs->run_threads = 0;
            } else  { // num_bytes == 0
                timeouts += TIMEOUT_USEC;
                if(timeouts >= 60*1000*1000) {
                    fprintf(stdout, "No packets received for 60 seconds on port %d.\n", bs->port);
                    syslog(LOG_WARNING, "no packets received for 60 seconds on port %d\n", bs->port);
                    // Send self the INT signal (simulate ctrl-c)
                    raise(SIGINT);
                    timeouts = 0;
                }
            }
            //printf("run_threads=%d num_bytes=%d errno=%d(%d)\n", bs->run_threads, num_bytes, errno, EINTR);
        }

        if (num_bytes > 0) {
            //printf("Wrote in a packet: size=%d slot=%d\n", num_bytes, this_slot - bs->buf->list_ptr);
            this_slot->size = num_bytes;
            // Mark this slot ready for readout
            sem_post(&this_slot->read_mutex);
            bs->buf->write_ptr = this_slot->next;
        }
    }
    close(sock);
    return NULL;
}

socket_t setup_network_listener(short port) {
    /* Open up a UDP socket on the specified port for receiving data */
    int sock = -1;
    struct sockaddr_in my_addr; // server's address information

    // create a new UDP socket descriptor
    sock = socket(PF_INET, SOCK_DGRAM, 0);
    if (sock == -1) {
        perror(__FILE__ " socket");
        return -1;
    }

    // initialize local address struct
    my_addr.sin_family = AF_INET; // host byte order
    my_addr.sin_port = htons(port); // short, network byte order
    my_addr.sin_addr.s_addr = htonl(INADDR_ANY); // listen on all interfaces
    memset(my_addr.sin_zero, 0, sizeof(my_addr.sin_zero));

    // bind socket to local address
    if (bind(sock, (SA *)&my_addr, sizeof(my_addr)) == -1){ 
        perror(__FILE__ " bind" );
        return -1;
    }   

    // prevent "address already in use" errors
    const int on = 1;
    if (setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, 
            (void *)&on, sizeof(on)) == -1) return -1;

    // Increase recv buffer for this sock
    int bufsize = 64*1024*1024; // 64 MB
    socklen_t ss = sizeof(int);
    int rv = setsockopt(sock, SOL_SOCKET, SO_RCVBUF, &bufsize, ss);
    if (rv<0) {
        perror("setsockopt");
    }
    rv = getsockopt(sock, SOL_SOCKET, SO_RCVBUF, &bufsize, &ss);
    if (rv==0) {
        printf("Using SO_RCVBUF size %d\n", bufsize);
    } else {
        perror("getsockopt");
    }

    return sock;
}
