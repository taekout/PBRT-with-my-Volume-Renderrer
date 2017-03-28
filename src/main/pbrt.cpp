// main/pbrt.cpp*
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"

// main program
int main(int argc, char *argv[]) {
    Options options;
    vector<string> filenames;
	string filename = "";
	bool bNonQuickRender = false;
    // Process command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--ncores")) options.nCores = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--outfile")) options.imageFile = argv[++i];
        else if (!strcmp(argv[i], "--quick")) options.quickRender = true;
		else if (!strcmp(argv[i], "--doublequick")) { options.doubleQuickRender = true; options.quickRender = false; }
		else if (!strcmp(argv[i], "--nonquick")) { bNonQuickRender = true; options.doubleQuickRender = false; options.quickRender = false; }
        else if (!strcmp(argv[i], "--quiet")) options.quiet = true;
        else if (!strcmp(argv[i], "--verbose")) options.verbose = true;
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            printf("usage: pbrt [--ncores n] [--outfile filename] [--quick] [--quiet] "
                   "[--verbose] [--help] <filename.pbrt> ...\n");
            return 0;
        }
        else {
			filenames.push_back(argv[i]);
			filename = argv[i];
		}
    }
	if(!bNonQuickRender && options.quickRender == false && options.doubleQuickRender == false) {
		options.quickRender = true;
		options.doubleQuickRender = false;
		printf("Quick Render by default.\n");
	}
	
    // Print welcome banner
    if (!options.quiet) {
        printf("pbrt version %s of %s at %s [Detected %d core(s)]\n",
               PBRT_VERSION, __DATE__, __TIME__, NumSystemCores());
        printf("Copyright (c)1998-2010 Matt Pharr and Greg Humphreys.\n");
        printf("The source code to pbrt (but *not* the book contents) is covered by the GNU GPL.\n");
        printf("See the file COPYING.txt for the conditions of the license.\n");
        fflush(stdout);
    }
    pbrtInit(options);
    // Process scene description
    PBRT_STARTED_PARSING();
    if (filenames.size() == 0) {
        // Parse scene from standard input
        ParseFile("-");
    } else {
        // Parse scene from input files
        for (u_int i = 0; i < filenames.size(); i++)
            if (!ParseFile(filenames[i]))
                Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
    }
    pbrtCleanup();
	/*if(filename.empty() == false) {
		string command;
		command = "\"C:\\Users\\shin\\Desktop\\Research Codes\\pbrt-v2\\bin\\exrtotiff.exe\"  ";
		command += "\"C:\\Users\\shin\\Desktop\\Research Codes\\pbrt-v2\\src\\pbrt.vs2010\\"; 
		command += filename.substr(13, filename.length() - 5 - 12); command += "exr\" ";
		command += "\"C:\\Users\\shin\\Desktop\\Research Codes\\pbrt-v2\\bin\\"; 
		command += filename.substr(13, filename.length() - 5 - 12); command += "tiff\"";
		printf("%s", command.c_str());
		//system(command.c_str());
	}*/
    return 0;
}


