#include <R.h> /* required for R specific stuff */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("qmgcv", String)
#else
#define _(String) (String)
#endif




