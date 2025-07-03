/** Get from https://gist.github.com/ocxtal/b2b4fab0373cae118bf98da0aa9822d8 **/
#ifndef _BITCNT_H
#define _BITCNT_H

#include <stdint.h>

/* Architecture-specific includes */
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
  #include <x86intrin.h>
  #define ARCH_X86
#elif defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
  /* ARM64/AArch64 - use GCC builtins */
  #define ARCH_ARM64
#endif

/**
 * misc bit operations (popcnt, tzcnt, and lzcnt)
 */

/**
 * @macro popcnt
 */
#if defined(ARCH_ARM64)
	/* ARM64: Use GCC builtin */
	#define popcnt(x)		( (uint64_t)__builtin_popcountll(x) )
#elif defined(__POPCNT__) && defined(ARCH_X86)
	/* x86 with POPCNT instruction */
	#define popcnt(x)		( (uint64_t)_mm_popcnt_u64(x) )
#else
	/* Software implementation for compatibility */
	// #warning "popcnt instruction is not enabled, using software implementation."
	static inline
	uint64_t popcnt(uint64_t n)
	{
		uint64_t c = 0;
		c = (n & 0x5555555555555555) + ((n>>1) & 0x5555555555555555);
		c = (c & 0x3333333333333333) + ((c>>2) & 0x3333333333333333);
		c = (c & 0x0f0f0f0f0f0f0f0f) + ((c>>4) & 0x0f0f0f0f0f0f0f0f);
		c = (c & 0x00ff00ff00ff00ff) + ((c>>8) & 0x00ff00ff00ff00ff);
		c = (c & 0x0000ffff0000ffff) + ((c>>16) & 0x0000ffff0000ffff);
		c = (c & 0x00000000ffffffff) + ((c>>32) & 0x00000000ffffffff);
		return(c);
	}
#endif

/**
 * @macro tzcnt
 * @brief trailing zero count (count #continuous zeros from LSb)
 */
#if defined(ARCH_ARM64)
	/* ARM64: Use GCC builtin */
	static inline
	uint64_t tzcnt(uint64_t n)
	{
		return (n == 0) ? 64 : (uint64_t)__builtin_ctzll(n);
	}
#elif defined(__BMI__) && defined(ARCH_X86)
	/* x86 with BMI instruction set */
	#define tzcnt(x)		( (uint64_t)_tzcnt_u64(x) )
#else
	/* Software implementation */
	// #warning "tzcnt instruction is not enabled, using software implementation."
	static inline
	uint64_t tzcnt(uint64_t n)
	{
		#if defined(__POPCNT__) && defined(ARCH_X86)
			return(popcnt(~n & (n - 1)));
		#elif defined(ARCH_X86)
			if(n == 0) {
				return(64);
			} else {
				int64_t res;
				__asm__( "bsfq %1, %0" : "=r"(res) : "r"(n) );
				return(res);
			}
		#else
			/* Generic bit manipulation implementation */
			if(n == 0) return 64;
			uint64_t count = 0;
			while((n & 1) == 0) {
				n >>= 1;
				count++;
			}
			return count;
		#endif
	}
#endif

/**
 * @macro lzcnt
 * @brief leading zero count (count #continuous zeros from MSb)
 */
#if defined(ARCH_ARM64)
	/* ARM64: Use GCC builtin */
	static inline
	uint64_t lzcnt(uint64_t n)
	{
		return (n == 0) ? 64 : (uint64_t)__builtin_clzll(n);
	}
#elif defined(__LZCNT__) && defined(ARCH_X86)
	/* x86 with LZCNT instruction */
	#define lzcnt(x)		( (uint64_t)_lzcnt_u64(x) )
#else
	/* Software implementation */
	// #warning "lzcnt instruction is not enabled, using software implementation."
	static inline
	uint64_t lzcnt(uint64_t n)
	{
		#if defined(ARCH_X86)
			if(n == 0) {
				return(64);
			} else {
				int64_t res;
				__asm__( "bsrq %1, %0" : "=r"(res) : "r"(n) );
				return(63 - res);
			}
		#else
			/* Generic bit manipulation implementation */
			if(n == 0) return 64;
			uint64_t count = 0;
			uint64_t mask = 0x8000000000000000ULL;
			while((n & mask) == 0) {
				mask >>= 1;
				count++;
			}
			return count;
		#endif
	}
#endif

#endif /* #ifndef _BITCNT_H */
