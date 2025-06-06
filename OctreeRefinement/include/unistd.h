/* ========================================================================== */
/*                                                                            */
/*   unistd.h (patched)                                                       */
/*   (c) 2012 Author                                                          */
/*                                                                            */
/*   Description:                                                             */
/*   Drop-in replacement for unistd.h on Windows; on other platforms,         */
/*   forward to the system <unistd.h>.                                        */
/*                                                                            */
/* ========================================================================== */

#ifndef _UNISTD_H
#define _UNISTD_H    1

#ifdef _WIN32

  /* This block is only used when compiling on Windows */
  #include <stdlib.h>
  #include <io.h>
  #include <getopt.h>   /* getopt at: https://gist.github.com/ashelly/7776712 */
  #include <process.h>  /* for getpid() and the exec..() family */
  #include <direct.h>   /* for _getcwd() and _chdir() */

  #define srandom srand
  #define random  rand

  /* Values for the second argument to access.
     These may be OR’d together.  */
  #define R_OK    4       /* Test for read permission.  */
  #define W_OK    2       /* Test for write permission.  */
  //#define X_OK  1       /* execute permission – unsupported in Windows */
  #define F_OK    0       /* Test for existence.  */

  #define access   _access
  #define dup2     _dup2
  #define execve   _execve
  #define ftruncate _chsize
  #define unlink   _unlink
  #define fileno   _fileno
  #define getcwd   _getcwd
  #define chdir    _chdir
  #define isatty   _isatty
  #define lseek    _lseek

  /* read, write, and close are NOT being #defined here, because while
     there are file‐handle‐specific versions for Windows, they probably
     don’t work for sockets.  You need to look at your app and consider
     whether to call e.g. closesocket(). */

  /* ssize_t is not defined on Windows by default: */
  #define ssize_t int

  #define STDIN_FILENO  0
  #define STDOUT_FILENO 1
  #define STDERR_FILENO 2

  /* Should be in some equivalent to <sys/types.h> */
  typedef __int16           int16_t; 
  typedef __int32           int32_t;
  typedef __int64           int64_t;
  typedef unsigned __int8   uint8_t;
  typedef unsigned __int16  uint16_t;
  typedef unsigned __int32  uint32_t;
  typedef unsigned __int64  uint64_t;

#else

  /* On Unix‐like systems (macOS, Linux, etc.), simply include the system’s <unistd.h> */
  #include <unistd.h>

#endif /* _WIN32 */

#endif /* _UNISTD_H */