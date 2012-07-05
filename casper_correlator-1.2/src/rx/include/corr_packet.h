#ifndef CORR_PACKET_H
#define CORR_PACKET_H

#include <netinet/in.h>
#include "packet_error.h"

#ifndef htonll
#ifdef _BIG_ENDIAN
#define htonll(x)   ((uint64_t)x)
#define ntohll(x)   ((uint64_t)x)
#else
#define htonll(x)   ((((uint64_t)htonl(x)) << 32) + htonl(((uint64_t)x) >> 32))
#define ntohll(x)   ((((uint64_t)ntohl(x)) << 32) + ntohl(((uint64_t)x) >> 32))
#endif
#endif

#define MAX_HEADER_SIZE     1024
#define MAX_PAYLOAD_SIZE    8192
#define OPTION_SIZE         8

#define PKT_ID              19282       // 0x4B52
#define TIME_PKT_ID         21587
#define PKT_VERSION         3
#define PKT_NOPTIONS        6

#define OPTID_SHIFT 	    48
#define OPTID_INSTIDS	    50
#define OPTID_TIMESTAMP     3
#define OPTID_HEAPLEN       4
#define OPTID_HEAPOFF       5
#define OPTID_CURRERR       52
#define OPTID_PKTINFO	    51
#define OPTID_HEAPPOINTER   53

#define PKT_DTYPE_UNSIGNED  0
#define PKT_DTYPE_SIGNED    1
#define PKT_DTYPE_FLOAT     2

// Define any options that need special packing/unpacking
/*typedef struct {
    uint16_t inst_id;
    int is_complex;
    int msb_first;
    uint8_t data_type;
    uint8_t bits_per_val;
    uint16_t vals_per_item;
} InstType;
*/

typedef struct {
 uint16_t instids_id;
 uint16_t instrument_id;
 uint16_t instance_id;
 uint16_t engine_id;
} InstIds;

// instrument_id of 3 is used for PAPER_FPGA X Engine
// instrument_id of 4 is used for PAPER_GPU X Engine
#define INSTRUMENT_ID_PAPER_FPGA_X_ENGINE (3)
#define INSTRUMENT_ID_PAPER_GPU_X_ENGINE (4)

//void unpack_insttype(InstType *it, uint64_t opt);
//uint64_t pack_insttype(InstType it);

typedef struct {
    uint32_t instance_id;
    uint32_t engine_id;
} InstSrc;

typedef struct {
    uint32_t len;
    uint32_t count;
} PktInfo;

void unpack_instsrc(InstSrc *is, uint64_t opt);
uint64_t pack_instsrc(InstSrc is);

// Define the packet
typedef struct {
    uint16_t n_options;         // Number of options to follow
    // Begin options = 16b identifier + 48b data
    //InstType insttype;          // Instrument Type
    //InstSrc instsrc;            // Instance ID + Engine #
    InstIds  instids;		// instrument, instance and engine IDs
    PktInfo  pktinfo;
    uint64_t timestamp;         // Timestamp
    uint64_t heap_pointer;          // NULL
    uint64_t heap_len;           // Data payload length
    uint64_t heap_off;           // Data payload offset
    uint64_t currerr;            // an error has occurred in the current packet
    uint64_t sync_time;         // the time in ADC samples since Unix Epoch  at which the correlator was last synced
    // And the payload
    char data[MAX_PAYLOAD_SIZE];
} CorrPacket;

#define HEADER_SIZE(pkt) ((pkt.n_options + 1) * OPTION_SIZE)

void init_packet(CorrPacket *pkt);
void unpack_header(CorrPacket *pkt, char *data);
void unpack_data(CorrPacket *pkt, char *data);
void unpack(CorrPacket *pkt, char *data);
void pack_header(CorrPacket pkt, char *data);
void pack_data(CorrPacket pkt, char *data);
void pack(CorrPacket pkt, char *data);
void unpack_instids(InstIds *id, uint64_t opt);
#endif
