#include "platform.h"
#include "internals.h"

posit64_t i32_to_p64(int32_t iA) {
    return convertDoubleToP64((double)iA);
}
