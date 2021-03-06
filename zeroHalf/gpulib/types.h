/**
 * @file   types.h
 *
 * @brief  Aliases for common data types
 *
 * @author Eyder Rios
 * @date   2011-09-12
 */

#ifndef __types_h
#define __types_h

#pragma GCC diagnostic ignored "-Wunused-result"

#ifdef __cplusplus

 #define EXTERN_C         extern "C"
 #define EXTERN_C_BEGIN   extern "C" {
 #define EXTERN_C_END     }

#else

 #define EXTERN_C
 #define EXTERN_C_BEGIN
 #define EXTERN_C_END

#endif


typedef unsigned char       byte,   *pbyte;     ///< Alias for unsigned char
typedef unsigned short      ushort, *pushort;   ///< Alias for unsigned short
typedef unsigned int        uint,   *puint;     ///< Alias for unsigned int
typedef unsigned long       ulong,  *pulong;    ///< Alias for unsigned long
typedef long long int       llong,  *pllong;    ///< Alias for long long
typedef unsigned long long  ullong, *pullong;   ///< Alias for unsigned long long

typedef void               *pointer;		    ///< Generic pointer (void *)

/*!
 * variant
 *
 * A variant data type.
 */
typedef union {
    char    v_char;
    byte    v_byte;
    short   v_short;
    ushort  v_ushort;
    int     v_int;
    uint    v_uint;
    long    v_long;
    ulong   v_ulong;
    llong   v_llong;
    ullong  v_ullong;
} variant;

/*!
 * vpointer
 *
 * A variant pointer.
 */
typedef union {
    const
    void        *p_cvoid;
    void        *p_void;
    char        *p_char;
    byte        *p_byte;
    short       *p_short;
    ushort      *p_ushort;
    int         *p_int;
    uint        *p_uint;
    long        *p_long;
    ulong       *p_ulong;
    llong       *p_llong;
    ullong      *p_ullong;
} vpointer;

/*!
 * To convert symbol to string
 */
#define SYM2STR(tok)        #tok
#define VAL2STR(tok)        SYM2STR(tok)

#endif
