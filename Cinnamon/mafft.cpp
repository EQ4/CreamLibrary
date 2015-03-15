#include "../c.library.h"

#include "OptimFFT.h"

static t_class *demifft_tilde_class;

typedef struct _demifft_tilde
{
    t_object            obj;
    t_sample            f;

    int                 partSize;
    t_sample*           buf;

    OptimFFT*           FFT;
    OptimFFT*           demiFFT;

    t_outlet*           out1_1;
    t_outlet*           out1_2;
    t_outlet*           out2_1;
    t_outlet*           out2_2;
    t_outlet*           out3_1;
    t_outlet*           out3_2;
    t_outlet*           out4_1;
    t_outlet*           out4_2;
    t_outlet*           out5_1;
    t_outlet*           out5_2;
}t_demifft_tilde;

static void *demifft_tilde_new()
{
	t_demifft_tilde *x   = (t_demifft_tilde *)pd_new(demifft_tilde_class);

    int   vectorSize = sys_getblksize();
    int huitiemeSize = vectorSize>>3;
    int    quartSize = vectorSize>>2;

    x->partSize = quartSize;

    x->demiFFT  = new OptimFFT(huitiemeSize);
    x->    FFT  = new OptimFFT(quartSize);

    x->buf      = new float[vectorSize];

    x->out1_1   = outlet_new(&x->obj, &s_signal);
    x->out1_2   = outlet_new(&x->obj, &s_signal);
    x->out2_1   = outlet_new(&x->obj, &s_signal);
    x->out2_2   = outlet_new(&x->obj, &s_signal);
    x->out3_1   = outlet_new(&x->obj, &s_signal);
    x->out3_2   = outlet_new(&x->obj, &s_signal);
    x->out4_1   = outlet_new(&x->obj, &s_signal);
    x->out4_2   = outlet_new(&x->obj, &s_signal);
    x->out5_1   = outlet_new(&x->obj, &s_signal);
    x->out5_2   = outlet_new(&x->obj, &s_signal);

    return (void *)x;
}

static void demifft_tilde_free(t_demifft_tilde *x)
{
    if (x->FFT)
    {
        delete x->FFT;
    }
    if (x->demiFFT)
    {
        delete x->demiFFT;
    }
    if (x->buf)
    {
        delete[] x->buf;
    }

    x->partSize   = 0;

	outlet_free(x->out1_1);
	outlet_free(x->out1_2);
	outlet_free(x->out2_1);
	outlet_free(x->out2_2);
	outlet_free(x->out3_1);
	outlet_free(x->out3_2);
	outlet_free(x->out4_1);
	outlet_free(x->out4_2);
	outlet_free(x->out5_1);
	outlet_free(x->out5_2);
}

static t_int *demifft_tilde_perform(t_int *w)
{
    t_demifft_tilde* x      = (t_demifft_tilde *)(w[1 ]);
    t_sample*        in     =        (t_sample *)(w[2 ]);
    t_sample*        out1_1 =        (t_sample *)(w[3 ]);
    t_sample*        out1_2 =        (t_sample *)(w[4 ]);
    t_sample*        out2_1 =        (t_sample *)(w[5 ]);
    t_sample*        out2_2 =        (t_sample *)(w[6 ]);
    t_sample*        out3_1 =        (t_sample *)(w[7 ]);
    t_sample*        out3_2 =        (t_sample *)(w[8 ]);
    t_sample*        out4_1 =        (t_sample *)(w[9 ]);
    t_sample*        out4_2 =        (t_sample *)(w[10]);
    t_sample*        out5_1 =        (t_sample *)(w[11]);
    t_sample*        out5_2 =        (t_sample *)(w[12]);
    int              n      =               (int)(w[13]);

    ::memcpy(x->buf, in, sys_getblksize()*sizeof(t_sample));

    ::memset(out1_1, 0, n*sizeof(t_sample)); ::memset(out1_2, 0, n*sizeof(t_sample));
    ::memset(out2_1, 0, n*sizeof(t_sample)); ::memset(out2_2, 0, n*sizeof(t_sample));
    ::memset(out3_1, 0, n*sizeof(t_sample)); ::memset(out3_2, 0, n*sizeof(t_sample));
    ::memset(out4_1, 0, n*sizeof(t_sample)); ::memset(out4_2, 0, n*sizeof(t_sample));
    ::memset(out5_1, 0, n*sizeof(t_sample)); ::memset(out5_2, 0, n*sizeof(t_sample));

    //FFT normale de 16 échantillons
    x->demiFFT->    FFT(out1_1, out1_2, x->buf            );
    //FFT normale de 16 échantillons (16 échantillons plus loin)
    x->demiFFT->    FFT(out2_1, out2_2, x->buf+x->partSize);
    //FFT normale de 32 échantillons
    x->    FFT->    FFT(out3_1, out3_2, x->buf            );
    //traitement qui va extraire les 2 premières FFTs et va les recombiner pour avoir les bins pairs
    x->demiFFT->EvenFFT(out4_1, out4_2, x->buf            );
    //traitement qui va extraire les 2 premières FFTs et va les recombiner pour avoir les bins impairs
    x->demiFFT->OddFFT(out5_1, out5_2, x->buf            );

    return (w+14);
}

static void demifft_tilde_dsp(t_demifft_tilde *x, t_signal **sp)
{
    dsp_add(demifft_tilde_perform, 13, x, sp[0]->s_vec, sp[1]->s_vec, sp[2 ]->s_vec, sp[3]->s_vec, sp[4]->s_vec,
                                                        sp[5]->s_vec, sp[6 ]->s_vec, sp[7]->s_vec, sp[8]->s_vec,
                                                        sp[9]->s_vec, sp[10]->s_vec, sp[0]->s_n);
}

extern "C" void demifft_tilde_setup(void) {
    demifft_tilde_class = class_new(gensym("demifft~"), (t_newmethod)demifft_tilde_new,
        (t_method)demifft_tilde_free, sizeof(t_demifft_tilde), CLASS_DEFAULT, (t_atomtype)0);

    class_addmethod(demifft_tilde_class, (t_method)demifft_tilde_dsp, gensym("dsp"), (t_atomtype)0);
    CLASS_MAINSIGNALIN(demifft_tilde_class, t_demifft_tilde, f);
}
