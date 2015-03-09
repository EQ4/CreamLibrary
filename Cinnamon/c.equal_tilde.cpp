/*
 * CicmWrapper
 *
 * A wrapper for Pure Data
 *
 * Copyright (C) 2013 Pierre Guillot, CICM - Université Paris 8
 * All rights reserved.
 *
 * Website  : http://www.mshparisnord.fr/HoaLibrary/
 * Contacts : cicm.mshparisnord@gmail.com
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Library General Public License as published
 * by the Free Software Foundation; either version 2 of the License.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

#include "../c.library.h"

typedef struct  _equal_tilde
{
    t_edspobj   j_obj;
    t_sample    f_value;
} t_equal_tilde;

t_eclass *equal_tilde_class;

extern void *equal_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_equal_tilde *x = (t_equal_tilde *)eobj_new(equal_tilde_class);
    x->f_value       = 0;
    eobj_dspsetup((t_ebox *)x, 2, 1);
    return (x);
}

extern void equal_tilde_perform_vector(t_equal_tilde *x, t_object *d, t_sample **ins, long ni, t_sample **outs, long no, long sampleframes, long f,void *up)
{
    while(--sampleframes)
    {
        outs[0][sampleframes] = t_sample(ins[0][sampleframes] == ins[1][sampleframes]);
    }
}

extern void equal_tilde_perform_scalar(t_equal_tilde *x, t_object *d, t_sample **ins, long ni, t_sample **outs, long no, long sampleframes, long f,void *up)
{
    while(--sampleframes)
    {
        outs[0][sampleframes] = t_sample(ins[0][sampleframes] == x->f_value);
    }
}

extern void equal_tilde_dsp(t_equal_tilde *x, t_object *dsp, short *count, double samplerate, long maxvectorsize, long flags)
{
    if(count[1])
    {
        object_method(dsp, gensym("dsp_add"), x, (method)equal_tilde_perform_vector, 0, NULL);
    }
    else
    {
        object_method(dsp, gensym("dsp_add"), x, (method)equal_tilde_perform_scalar, 0, NULL);
    }
}

extern t_pd_err equal_tilde_notify(t_equal_tilde *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
    if(msg == gensym("attr_modified"))
    {
        if(s == gensym("value"))
        {
            post("Your value has changed : %f", x->f_value);
        }
    }
    return 0;
}

extern "C"  void setup_c0x2eequal_tilde(void)
{
    t_eclass *c;
    c = eclass_new("c.==~", (method)equal_tilde_new, (method)eobj_dspfree, (short)sizeof(t_equal_tilde), 0L, A_GIMME, 0);
    
    eclass_dspinit(c);
    cream_initclass(c);
    
    eclass_addmethod(c, (method)equal_tilde_dsp,        "dsp",              A_NULL, 0);
    eclass_addmethod(c, (method)equal_tilde_notify,     "notify",           A_NULL, 0);
    eclass_register(CLASS_OBJ, c);
    
    CLASS_ATTR_FLOAT                (c, "value", 0, t_equal_tilde, f_value);
    CLASS_ATTR_LABEL                (c, "value", 0, "Comparative Value");
    CLASS_ATTR_ORDER                (c, "value", 0, "1");
    CLASS_ATTR_DEFAULT              (c, "value", 0, "0.");
    CLASS_ATTR_STYLE                (c, "value", 0, "number");
    
    equal_tilde_class = c;
}







