#ifndef SDISP_H_
#define SDISP_H_

#define SDISP_MAGIC 0x5344

#define SDISP_POWER_SPECTRUM 0x0001
#define SDISP_PHASE_SPECTRUM 0x0002
#define MAX_SDISP_IP_COUNT 10
#define SDISP_PACKET_SIZE 8192

struct sdisp_complex {
 float re;
 float im;
} __attribute__ ((packed));

struct sdisp_header{
    uint16_t magic_number;
    uint16_t corr_product_id;
    uint16_t offset;
    uint16_t length;
    uint64_t timestamp_ms;
} __attribute__ ((packed));

#endif
