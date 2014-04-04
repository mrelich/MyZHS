#include <features.h>
#if __GLIBC_MINOR__ > 0

#include <fenv.h>

#else

#include <fpu_control.h>

#endif

static void __attribute__ ((constructor))
     trapfpe ()
{
#if __GLIBC_MINOR__ > 0

  feenabletraps ((FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));

#else

  __setfpucw (_FPU_DEFAULT &
	      ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM));

#endif
}

void fpenab_ () {
  trapfpe();
  }

/*
  Enable exception traps for the exceptions specified by traps.
  */
int feenabletraps(int traps)
{
  unsigned short cw;

  asm volatile ("fnstcw %0":"=m" (cw));
  cw &= ~traps;
  asm volatile ("fldcw %0"::"m" (cw));
  return 0;
}









