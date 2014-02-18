#ifndef VERSION_H
#define VERSION_H

namespace AutoVersion{

	//Date Version Types
	static const char DATE[] = "03";
	static const char MONTH[] = "07";
	static const char YEAR[] = "2008";
	static const double UBUNTU_VERSION_STYLE = 8.07;

	//Software Status
	static const char STATUS[] = "Alpha";
	static const char STATUS_SHORT[] = "a";

	//Standard Version Type
	static const long MAJOR = 0;
	static const long MINOR = 0;
	static const long BUILD = 0;
	static const long REVISION = 0;

	//Miscellaneous Version Types
	static const long BUILDS_COUNT = 265;
	#define RC_FILEVERSION 0,0,0,0
	#define RC_FILEVERSION_STRING "0, 0, 0, 0\0"
	static const char FULLVERSION_STRING[] = "0.0.0.0";

	//These values are to keep track of your versioning state, don't modify them.
	static const long BUILD_HISTORY = 0;

    // #define VERBOSE_FLAG
    // #define VERBOSE_OUT
    #define GNUPLOT_SUPPORT
    #define LATEX_SUPPORT
}
#endif //VERSION_h
