/** Get from https://gist.github.com/ocxtal/b2b4fab0373cae118bf98da0aa9822d8 **/
#ifndef _BITCNT_H
#define _BITCNT_H

#include <x86intrin.h>
#include <stdint.h>

/**
 * misc bit operations (popcnt, tzcnt, and lzcnt)
 */

/**
 * @macro popcnt
 */
#ifdef __POPCNT__
	#define popcnt(x)		( (uint64_t)_mm_popcnt_u64(x) )
#else
	// #warning "popcnt instruction is not enabled."
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
#ifdef __BMI__
	/** immintrin.h is already included */
	#define tzcnt(x)		( (uint64_t)_tzcnt_u64(x) )
#else
	// #warning "tzcnt instruction is not enabled."
	static inline
	uint64_t tzcnt(uint64_t n)
	{
		#ifdef __POPCNT__
			return(popcnt(~n & (n - 1)));
		#else
			if(n == 0) {
				return(64);
			} else {
				int64_t res;
				__asm__( "bsfq %1, %0" : "=r"(res) : "r"(n) );
				return(res);
			}
		#endif
	}
#endif

/**
 * @macro lzcnt
 * @brief leading zero count (count #continuous zeros from MSb)
 */
#ifdef __LZCNT__
	#define lzcnt(x)		( (uint64_t)_lzcnt_u64(x) )
#else
	// #warning "lzcnt instruction is not enabled."
	static inline
	uint64_t lzcnt(uint64_t n)
	{
		if(n == 0) {
			return(64);
		} else {
			int64_t res;
			__asm__( "bsrq %1, %0" : "=r"(res) : "r"(n) );
			return(63 - res);
		}
	}
#endif

#endif /* #ifndef _BITCNT_H */
