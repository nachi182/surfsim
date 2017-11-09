#ifndef MYCONST_H
#define MYCONST_H

// Allowed leakage reduction techniques
enum corTech {no, quick, gate, circuit};

// Can leakage be detected? (modifies the decoding graph)
static const bool lkgDetect = true;

// Set to true to print debugging messages on the screen
static const bool debug = false;

#endif
