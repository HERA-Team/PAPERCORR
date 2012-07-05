#include <stdio.h>
#include "include/corr_packet.h"

void unpack_instids(InstIds *id, uint64_t opt) {
 id->instrument_id = (uint16_t)((opt >> 32) & 0xFFFF);
 id->instance_id = (uint16_t)((opt >> 16) & 0xFFFF);
 id->engine_id = (uint16_t)((opt >> 0) & 0xFFFF);
}

uint64_t pack_instids(InstIds id) {
 return (((uint64_t)id.instrument_id << 32) +
	 ((uint64_t)id.instance_id << 16) + 
	 ((uint64_t)id.engine_id << 0));
}

void unpack_pktinfo(PktInfo *pi, uint64_t opt) {
    pi->len = (uint32_t)((opt >> 24) & 0xFFFFFF);
    pi->count = (uint32_t)((opt >> 0) & 0xFFFFFF);
}

uint64_t pack_pktinfo(PktInfo pi) {
 return (((uint64_t)pi.len << 24) + ((uint64_t)pi.count << 0));
}

void init_packet(CorrPacket *pkt) {
    pkt->n_options = 6;
    pkt->instids.instrument_id = 3;
    pkt->instids.instance_id = 0;
    pkt->instids.engine_id = 0;
    pkt->pktinfo.len = 0;
    pkt->pktinfo.count = 0;
    pkt->timestamp = 0;
    pkt->heap_len = 0;
    pkt->heap_off = 0;
    pkt->currerr = 0;
    pkt->heap_pointer = 0;
    pkt->sync_time = 0;
}

#define OPT(data,n) ntohll(((uint64_t *)(data + n * OPTION_SIZE))[0])
#define OPTID(opt) ((uint16_t) 0xFFFF & (opt >> OPTID_SHIFT))
#define OPTVAL(opt) (opt & 0xFFFFFFFFFFFF)

void unpack_header(CorrPacket *pkt, char *data) {
    // Unpack header from data buffer.  Assumes data is a valid pointer.
    uint64_t opt;
    time_t ppt;
    opt = OPT(data,0);  // Endian conversion happens here
    // Make sure this is our kind of packet
    if (((opt >> 48) & 0xFFFF) == TIME_PKT_ID) {
     ppt = OPT(data,1);
     //fprintf(stderr, "Got time sync packet: %lu, %s\n", OPT(data,1), ctime(&ppt));
     pkt->sync_time = OPT(data,1);
     pkt->pktinfo.len = 0;
     return;
    }

    if (((opt >> 48) & 0xFFFF) != PKT_ID) {
	    //fprintf(stderr, "Received ID %i\n", ((opt >> 48) & 0xFFFF));
        throw PacketError("Wrong PKT_ID... not a SPEAD packet.");
    }
    if (((opt >> 32) & 0xFFFF) != PKT_VERSION)
        throw PacketError("Wrong PKT_VERSION");
    pkt->n_options = (uint16_t) 0xFFFF & (opt >> 0);
    // Unpack options (and be agnostic about how many there are)
    int i;
    for (i=1; i <= pkt->n_options; i++) {
        opt = OPT(data,i);  // Endian conversion happens here
        switch(OPTID(opt)) {
            case OPTID_INSTIDS:
                unpack_instids(&(pkt->instids),OPTVAL(opt)); break;
            case OPTID_PKTINFO:
		        unpack_pktinfo(&(pkt->pktinfo), OPTVAL(opt)); break;
            case OPTID_TIMESTAMP:
                pkt->timestamp = OPTVAL(opt); break;
            case OPTID_HEAPOFF:
                pkt->heap_off = OPTVAL(opt); break;
            case OPTID_CURRERR:
                pkt->currerr = OPTVAL(opt); break;
	    case OPTID_HEAPPOINTER:
                pkt->heap_pointer = OPTVAL(opt); break;
            default: break;
        }
    }
    // Check the options we really can't do without
    if (pkt->pktinfo.len > MAX_PAYLOAD_SIZE)
        throw PacketError("pkt.pktinfo.len > MAX_PAYLOAD_SIZE");
}

#define IND8(data,n)  (((uint8_t  *)data)[n])
#define IND16(data,n) (((uint16_t *)data)[n])
#define IND32(data,n) (((int32_t  *)data)[n])
#define IND64(data,n) (((uint64_t *)data)[n])
#define FLOAT(data,n) (((float    *)data)[n])

#define GET8(data,n) IND8(data,n)
#define GET16(data,n) ntohs(IND16(data,n))
#define GET32(data,n) ntohl(IND32(data,n))
#define GET64(data,n) ntohll(IND64(data,n))

#define SWAP8(val) ((uint8_t)val)
#define SWAP16(val) ((((uint16_t)SWAP8(val)) << 8) + SWAP8(((uint16_t)val) >> 8))
#define SWAP32(val) ((((int32_t)SWAP16(val)) << 16) + SWAP16(((int32_t)val) >> 16))
#define SWAP64(val) ((((uint64_t)SWAP32(val)) << 32) + SWAP32(((uint64_t)val) >> 32))

void unpack_data(CorrPacket *pkt, char *data) {
    // Unpack payload from data buffer.
    // Assumes data pointer (pointing to start of payload) is valid, and 
    // packet header info has been set.  
    // Based solely on instrument type.
    uint64_t i;
    switch (pkt->instids.instrument_id) {
        case INSTRUMENT_ID_PAPER_FPGA_X_ENGINE:
            // 32 bit big endian (aka network byte order) integer data
            for (i=0; i < (pkt->pktinfo.len / 4); i++) {
                IND32(pkt->data,i) = GET32(data,i);
            }
            break;
        case INSTRUMENT_ID_PAPER_GPU_X_ENGINE:
            // 32 bit little endian float data
            for (i=0; i < (pkt->pktinfo.len / 4); i++) {
                FLOAT(pkt->data,i) = FLOAT(data,i);
            }
            break;
        default:
            char msg[128];
            sprintf(msg, "Unknown instrument id (%d)", pkt->instids.instrument_id);
            throw PacketError(msg);
    }
}

void unpack(CorrPacket *pkt, char *data) {
    unpack_header(pkt, data);
    unpack_data(pkt, data + HEADER_SIZE((*pkt)));
}

#define MKOPT(id,val) htonll(((uint64_t)id << OPTID_SHIFT) + val)

void pack_header(CorrPacket pkt, char *data) {
    // Pack information from CorrPacket header into data buffer
    IND64(data,0) = MKOPT(PKT_ID, (((uint64_t)PKT_VERSION << 32) + PKT_NOPTIONS));
    // Pack options
    IND64(data,1) = MKOPT(OPTID_INSTIDS, pack_instids(pkt.instids));
    IND64(data,2) = MKOPT(OPTID_TIMESTAMP, pkt.timestamp);
    IND64(data,3) = MKOPT(OPTID_CURRERR, pkt.currerr);
    IND64(data,4) = MKOPT(OPTID_PKTINFO, pack_pktinfo(pkt.pktinfo));
    IND64(data,5) = MKOPT(OPTID_HEAPOFF, pkt.heap_off);
    IND64(data,6) = MKOPT(OPTID_HEAPPOINTER, pkt.heap_pointer);
}

#define SET8(data,n) IND8(data,n)
#define SET16(data,n) htons(IND16(data,n))
#define SET32(data,n) htonl(IND32(data,n))
#define SET64(data,n) htonll(IND64(data,n))

void pack_data(CorrPacket pkt, char *data) {
    // Pack packet data into a data buffer
    uint64_t i;
    switch(pkt.instids.instrument_id) {
        case 3:
            for (i=0; i < (pkt.pktinfo.len / 4); i++) {
                IND32(data,i) = SET32(pkt.data,i);
            }
            break;
        default: break;
    }
}

void pack(CorrPacket pkt, char *data) {
    // Pack an entire packet into a data buffer
    pack_header(pkt, data);
    pack_data(pkt, data + HEADER_SIZE(pkt));
}

