/*
// Copyright (c) 2012-2014 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_CREAM_LIBRARY
#define DEF_CREAM_LIBRARY

extern "C"
{
#include "ThirdParty/PureData/Sources/cicm_wrapper.h"
}

void cream_initclass(t_eclass* c);

extern "C" void setup_c0x2elibrary(void);
extern "C" void libpd_loadcream(void);

// Cinnamon
extern "C"  void setup_c0x2eatan_tilde(void);
extern "C"  void setup_c0x2esin_tilde(void);
extern "C"  void setup_c0x2eequal_tilde(void);

// Caramel
extern "C" void setup_c0x2econvolve_tilde(void);
extern "C" void setup_c0x2ecount_tilde(void);
extern "C" void setup_c0x2efft_tilde(void);
extern "C" void setup_c0x2efreeverb_tilde(void);
extern "C" void setup_c0x2ematrix_tilde(void);

// Chocolate
extern "C" void setup_c0x2ebang(void);
extern "C" void setup_c0x2eblackboard(void);
extern "C" void setup_c0x2ebreakpoints(void);
extern "C" void setup_c0x2ecolorpanel(void);
extern "C" void setup_c0x2edsp_tilde(void);
extern "C" void setup_c0x2egain_tilde(void);
extern "C" void setup_c0x2eincdec(void);
extern "C" void setup_c0x2eknob(void);
extern "C" void setup_c0x2ematrix(void);
extern "C" void setup_c0x2emenu(void);
extern "C" void setup_c0x2emeter_tilde(void);
extern "C" void setup_c0x2enumber(void);
extern "C" void setup_c0x2enumber_tilde(void);
extern "C" void setup_c0x2eplane(void);
extern "C" void setup_c0x2epreset(void);
extern "C" void setup_c0x2eradio(void);
extern "C" void setup_c0x2erslider(void);
extern "C" void setup_c0x2escope_tilde(void);
extern "C" void setup_c0x2eslider(void);
extern "C" void setup_c0x2espectroscope_tilde(void);
extern "C" void setup_c0x2etab(void);
extern "C" void setup_c0x2etoggle(void);

// Coffee

extern "C" void setup_c0x2eloadmess(void);
extern "C" void setup_c0x2epak(void);
extern "C" void setup_c0x2epatcherargs(void);
extern "C" void setup_c0x2epatcherinfos(void);
extern "C" void setup_c0x2epatchermess(void);
extern "C" void setup_c0x2eprepend(void);
extern "C" void setup_c0x2emousestate(void);

#endif
