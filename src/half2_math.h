__device__ __forceinline__ __half2 operator+(const __half2 &lh, const __half2 &rh) { return __hadd2(lh, rh); }
__device__ __forceinline__ __half2 operator-(const __half2 &lh, const __half2 &rh) { return __hsub2(lh, rh); }
__device__ __forceinline__ __half2 operator/(const __half2 &lh, const __half2 &rh) { return __h2div(lh, rh); }
// __device__ __forceinline__ __half2 operator*(const __half2 &lh, const __half2 &rh) { return __hmul2(lh, rh); }
__device__ __inline__ __half2 operator*(const __half2 &lh, const __half2 &rh) { 
    return __half2{lh.x * rh.x - lh.y * rh.y, lh.x * rh.y + lh.y * rh.x};
}

__device__ __forceinline__ __half2& operator+=(__half2 &lh, const __half2 &rh) { lh = __hadd2(lh, rh); return lh; }
__device__ __forceinline__ __half2& operator-=(__half2 &lh, const __half2 &rh) { lh = __hsub2(lh, rh); return lh; }
__device__ __forceinline__ __half2& operator/=(__half2 &lh, const __half2 &rh) { lh = __h2div(lh, rh); return lh; }
// __device__ __forceinline__ __half2& operator*=(__half2 &lh, const __half2 &rh) { lh = __hmul2(lh, rh); return lh; }
__device__ __inline__ __half2& operator*=(__half2 &lh, const __half2 &rh) { lh = lh * rh; return lh; }

__device__ __forceinline__ __half2 &operator++(__half2 &h)      { __half2_raw one; one.x = 0x3C00U; one.y = 0x3C00U; h = __hadd2(h, one); return h; }
__device__ __forceinline__ __half2 &operator--(__half2 &h)      { __half2_raw one; one.x = 0x3C00U; one.y = 0x3C00U; h = __hsub2(h, one); return h; }
__device__ __forceinline__ __half2  operator++(__half2 &h, const int ignored)
{
    // ignored on purpose. Parameter only needed to distinguish the function declaration from other types of operators.
    static_cast<void>(ignored);

    const __half2 ret = h;
    __half2_raw one;
    one.x = 0x3C00U;
    one.y = 0x3C00U;
    h = __hadd2(h, one);
    return ret;
}
__device__ __forceinline__ __half2  operator--(__half2 &h, const int ignored)
{
    // ignored on purpose. Parameter only needed to distinguish the function declaration from other types of operators.
    static_cast<void>(ignored);

    const __half2 ret = h;
    __half2_raw one;
    one.x = 0x3C00U;
    one.y = 0x3C00U;
    h = __hsub2(h, one);
    return ret;
}

__device__ __forceinline__ __half2 operator+(const __half2 &h) { return h; }
__device__ __forceinline__ __half2 operator-(const __half2 &h) { return __hneg2(h); }

__device__ __forceinline__ bool operator==(const __half2 &lh, const __half2 &rh) { return __hbeq2(lh, rh); }
__device__ __forceinline__ bool operator!=(const __half2 &lh, const __half2 &rh) { return __hbneu2(lh, rh); }
__device__ __forceinline__ bool operator>(const __half2 &lh, const __half2 &rh) { return __hbgt2(lh, rh); }
__device__ __forceinline__ bool operator<(const __half2 &lh, const __half2 &rh) { return __hblt2(lh, rh); }
__device__ __forceinline__ bool operator>=(const __half2 &lh, const __half2 &rh) { return __hbge2(lh, rh); }
__device__ __forceinline__ bool operator<=(const __half2 &lh, const __half2 &rh) { return __hble2(lh, rh); }
