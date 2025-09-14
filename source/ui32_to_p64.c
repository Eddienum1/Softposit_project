#include <stdint.h>

#include "platform.h"
#include "internals.h"

posit64_t ui32_to_p64(uint32_t a) {
    return ui64_to_p64((uint64_t)a);
}
