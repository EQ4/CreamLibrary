/*
 * Cream Library
 * Copyright (C) 2013 Pierre Guillot, CICM - Université Paris 8
 * All rights reserved.
 * Website  : https://github.com/CICM/CreamLibrary
 * Contacts : cicm.mshparisnord@gmail.com
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 */

#include "../c.library.hpp"
#if (_MSC_VER >= 1800)
static double round(double val)
{
	return floor(val + 0.5);
}
#endif

typedef struct  _colorpanel
{
	t_ebox      j_box;
    t_outlet*   f_out_rgb;
    t_outlet*   f_out_hsl;
    t_outlet*   f_out_hex;
    t_hsla**    f_matrix_colorpanel;
    t_pt        f_matrix_sizes;
    char        f_reverse;
    float       f_saturation;
    float       f_hue;
    float       f_lightness;

    t_pt        f_color_picked;
    t_pt        f_color_hover;

	t_rgba		f_color_background;
	t_rgba		f_color_border;
} t_colorpanel;

static t_eclass *colorpanel_class;

static void *colorpanel_new(t_symbol *s, int argc, t_atom *argv)
{
    int i;
	t_colorpanel *x = (t_colorpanel *)eobj_new(colorpanel_class);
    t_binbuf* d = binbuf_via_atoms(argc,argv);

    if(x && d)
    {
        x->f_matrix_sizes.x = 1;
        x->f_matrix_sizes.y = 1;
        x->f_matrix_colorpanel   = (t_hsla **)malloc(1 * sizeof(t_hsla *));
        for(i = 0; i < 1; i++)
        {
            x->f_matrix_colorpanel[i]   = (t_hsla *)malloc(1 * sizeof(t_hsla));
        }
        
        ebox_new((t_ebox *)x, 0 | EBOX_GROWINDI);
        
        x->f_color_hover.x = -10;
        x->f_color_hover.y = -10;
        x->f_color_picked.x = -10;
        x->f_color_picked.y = -10;
        x->f_out_rgb = outlet_new((t_object *)x, &s_list);
        x->f_out_hsl = outlet_new((t_object *)x, &s_list);
        x->f_out_hex = outlet_new((t_object *)x, &s_symbol);
        
        ebox_attrprocess_viabinbuf(x, d);
        ebox_ready((t_ebox *)x);
        return (x);
    }
    
    return NULL;
}

static void colorpanel_free(t_colorpanel *x)
{
    int i;
    ebox_free((t_ebox *)x);
    for(i = 0; i < x->f_matrix_sizes.x; i++)
    {
        free(x->f_matrix_colorpanel[i]);
    }
    free(x->f_matrix_colorpanel);
}

static void colorpanel_output(t_colorpanel *x)
{
    t_rgba color_rgb;
    t_hsla color_hls;
    t_symbol* color_hex;
    t_atom av[3];
    if(x->f_color_picked.x >= 0 && x->f_color_picked.y >= 0)
    {
        t_pd* send = ebox_getsender((t_ebox *) x);
        color_hls = x->f_matrix_colorpanel[(int)x->f_color_picked.x][(int)x->f_color_picked.y];
        color_rgb = hsla_to_rgba(&color_hls);
        color_hex = gensym(rgba_to_hex(&color_rgb));
        atom_setfloat(av, color_rgb.red);
        atom_setfloat(av+1, color_rgb.green);
        atom_setfloat(av+2, color_rgb.blue);
        outlet_list(x->f_out_rgb, &s_list, 3, av);
        if(send)
        {
            pd_list(send, &s_list, 3, av);
        }
        atom_setfloat(av, color_hls.hue);
        atom_setfloat(av+1, color_hls.saturation);
        atom_setfloat(av+2, color_hls.lightness);
        outlet_list(x->f_out_hsl, &s_list, 3, av);
        outlet_symbol(x->f_out_hex, color_hex);
        
    }
}

static void colorpanel_set(t_colorpanel *x, t_symbol *s, int ac, t_atom *av)
{
    if(ac > 1 && av)
    {
        if(atom_gettype(av) == A_FLOAT)
            x->f_color_picked.x = pd_clip((int)atom_getfloat(av), 0, x->f_matrix_sizes.x-1);
        if(atom_gettype(av+1) == A_FLOAT)
            x->f_color_picked.y = pd_clip((int)atom_getfloat(av+1), 0, x->f_matrix_sizes.y-1);

        ebox_invalidate_layer((t_ebox *)x, cream_sym_picked_layer);
        ebox_redraw((t_ebox *)x);
    }
}

static void colorpanel_list(t_colorpanel *x, t_symbol *s, int ac, t_atom *av)
{
    if(ac > 1 && av)
    {
        if(atom_gettype(av) == A_FLOAT)
            x->f_color_picked.x = pd_clip((int)atom_getfloat(av), 0, x->f_matrix_sizes.x-1);
        if(atom_gettype(av+1) == A_FLOAT)
            x->f_color_picked.y = pd_clip((int)atom_getfloat(av+1), 0, x->f_matrix_sizes.y-1);

        colorpanel_output(x);
        ebox_invalidate_layer((t_ebox *)x, cream_sym_picked_layer);
        ebox_redraw((t_ebox *)x);
    }
}

static t_pd_err colorpanel_notify(t_colorpanel *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
	if(msg == cream_sym_attr_modified)
	{
        ebox_invalidate_layer((t_ebox *)x, cream_sym_background_layer);
        ebox_invalidate_layer((t_ebox *)x, cream_sym_hover_layer);
        ebox_invalidate_layer((t_ebox *)x, cream_sym_picked_layer);
	}
	return 0;
}

static void colorpanel_getdrawparams(t_colorpanel *x, t_object *patcherview, t_edrawparams *params)
{
    params->d_borderthickness   = 2;
    params->d_cornersize        = 2;
    params->d_bordercolor       = x->f_color_border;
    params->d_boxfillcolor      = x->f_color_background;
}

static void colorpanel_oksize(t_colorpanel *x, t_rect *newrect)
{
    float ratio;
    newrect->width = pd_clip_min(newrect->width, (float)x->f_matrix_sizes.x * 10.f);
    newrect->height = pd_clip_min(newrect->height, (float)x->f_matrix_sizes.y * 10.f);
    
    ratio = (newrect->width - 1.) / (float)x->f_matrix_sizes.x;
    if(ratio - (int)ratio != 0)
    {
        ratio = round(ratio);
        newrect->width = ratio * (float)x->f_matrix_sizes.x + 1.;
    }
    ratio = (newrect->height - 1.) / (float)x->f_matrix_sizes.y;
    if(ratio - (int)ratio != 0)
    {
        ratio = round(ratio);
        newrect->height = ratio * (float)x->f_matrix_sizes.y + 1.;
    }
    
    newrect->width = pd_clip_min(newrect->width, 30.);
    newrect->height = pd_clip_min(newrect->height, 10.);
}

static void draw_background(t_colorpanel *x, t_object *view, t_rect *rect)
{
	t_elayer *g = ebox_start_layer((t_ebox *)x, cream_sym_background_layer, rect->width, rect->height);
	if (g)
	{
        const float block_width = rect->width / (float)x->f_matrix_sizes.x;
        const float block_height = rect->height / (float)x->f_matrix_sizes.y;
        float incx = 0.f;
        for(int i = 0; i < x->f_matrix_sizes.x; i++)
        {
            float incY = 0.f;
            for(int j = 0; j < x->f_matrix_sizes.y; j++)
            {
                egraphics_set_color_hsla(g, &x->f_matrix_colorpanel[i][j]);
                egraphics_rectangle(g, incx + 1.f, incY + 1.f, block_width - 2.f, block_height - 2.f);
                egraphics_fill(g);
                incY += block_height;
            }
            incx += block_width;
        }

        ebox_end_layer((t_ebox*)x, cream_sym_background_layer);
	}
	ebox_paint_layer((t_ebox *)x, cream_sym_background_layer, 0, 0);
}

static void draw_picked(t_colorpanel *x, t_object *view, t_rect *rect)
{
	t_elayer *g = ebox_start_layer((t_ebox *)x, cream_sym_picked_layer, rect->width, rect->height);
	if(g)
	{
        const float block_width = rect->width / (float)x->f_matrix_sizes.x;
        const float block_height = rect->height / (float)x->f_matrix_sizes.y;
        if(x->f_color_picked.x >= 0 && x->f_color_picked.y >= 0)
        {
            egraphics_set_color_rgba(g, &x->f_color_border);
            egraphics_set_line_width(g, 2);
            egraphics_rectangle(g,  (float)x->f_color_picked.x * block_width + 1.f,
                                            (float)x->f_color_picked.y * block_height + 1.f,
                                            block_width - 2.f, block_height - 2.f);
            egraphics_stroke(g);
        }
        ebox_end_layer((t_ebox*)x, cream_sym_picked_layer);
	}
	ebox_paint_layer((t_ebox *)x, cream_sym_picked_layer, 0, 0);
}

static void draw_hover(t_colorpanel *x, t_object *view, t_rect *rect)
{
	t_elayer *g = ebox_start_layer((t_ebox *)x, cream_sym_hover_layer, rect->width, rect->height);
	if(g)
	{
        const float block_width = rect->width / (float)x->f_matrix_sizes.x;
        const float block_height = rect->height / (float)x->f_matrix_sizes.y;
        if(x->f_color_hover.x >= 0 && x->f_color_hover.y >= 0)
        {
            egraphics_set_color_hsla(g, &x->f_matrix_colorpanel[(int)x->f_color_hover.x][(int)x->f_color_hover.y]);
            egraphics_set_line_width(g, 2);
            egraphics_rectangle(g,  (float)x->f_color_hover.x * block_width,
                                    (float)x->f_color_hover.y * block_height,
                                    block_width, block_height);
            egraphics_fill(g);
        }
        
        ebox_end_layer((t_ebox*)x, cream_sym_hover_layer);
	}
	ebox_paint_layer((t_ebox *)x, cream_sym_hover_layer, 0, 0);
}

static void colorpanel_paint(t_colorpanel *x, t_object *view)
{
    t_rect rect;
    ebox_get_rect_for_view((t_ebox *)x, &rect);
    draw_background(x, view, &rect);
    draw_picked(x, view, &rect);
    draw_hover(x, view, &rect);
}

static void colorpanel_mousemove(t_colorpanel *x, t_object *patcherview, t_pt pt, long modifiers)
{
    x->f_color_hover.x = pd_clip((int)(pt.x / (x->j_box.b_rect.width / (float)x->f_matrix_sizes.x)), 0, x->f_matrix_sizes.x-1);
    x->f_color_hover.y = pd_clip((int)(pt.y / (x->j_box.b_rect.height / (float)x->f_matrix_sizes.y)), 0, x->f_matrix_sizes.y-1);
    ebox_invalidate_layer((t_ebox *)x, cream_sym_hover_layer);
    ebox_redraw((t_ebox *)x);
}

static void colorpanel_mousedown(t_colorpanel *x, t_object *patcherview, t_pt pt, long modifiers)
{
    x->f_color_hover.x = -10;
    x->f_color_hover.y = -10;
    x->f_color_picked.x = pd_clip((int)(pt.x / (x->j_box.b_rect.width / (float)x->f_matrix_sizes.x)), 0, x->f_matrix_sizes.x-1);
    x->f_color_picked.y = pd_clip((int)(pt.y / (x->j_box.b_rect.height / (float)x->f_matrix_sizes.y)), 0, x->f_matrix_sizes.y-1);
    ebox_invalidate_layer((t_ebox *)x, cream_sym_hover_layer);
    ebox_invalidate_layer((t_ebox *)x, cream_sym_picked_layer);
    ebox_redraw((t_ebox *)x);
    colorpanel_output(x);
}

static void colorpanel_mouseleave(t_colorpanel *x, t_object *patcherview, t_pt pt, long modifiers)
{
    x->f_color_hover.x = -10;
    x->f_color_hover.y = -10;
    ebox_invalidate_layer((t_ebox *)x, cream_sym_hover_layer);
    ebox_redraw((t_ebox *)x);
}

static void colorpanel_preset(t_colorpanel *x, t_binbuf *b)
{
   binbuf_addv(b, (char *)"sff", gensym("list"), x->f_color_picked.x, x->f_color_picked.y);
}

static void colorpanel_computecolors(t_colorpanel *x)
{
    int i, j;
    t_hsla color_ref;

    x->f_color_hover.x = -10;
    x->f_color_hover.y = -10;

    x->f_color_picked.x = -10;
    x->f_color_picked.y = -10;


    hsla_set(&color_ref, x->f_hue, x->f_saturation, x->f_lightness, 1.);
    if(x->f_reverse)
    {
        for(j = 0; j < x->f_matrix_sizes.y; j++)
        {
            color_ref.lightness = x->f_lightness;
            for (i = 0; i < x->f_matrix_sizes.x; i++)
            {
                x->f_matrix_colorpanel[i][j] = color_ref;
                color_ref.lightness -= (1. / (float)(x->f_matrix_sizes.x - 1));
                if(x->f_matrix_colorpanel[i][j].lightness > 1)
                    x->f_matrix_colorpanel[i][j].lightness = 1. - (x->f_matrix_colorpanel[i][j].lightness - 1.);
                else if(x->f_matrix_colorpanel[i][j].lightness < 0.)
                    x->f_matrix_colorpanel[i][j].lightness = -x->f_matrix_colorpanel[i][j].lightness;
            }
            color_ref.hue += (1. / (float)(x->f_matrix_sizes.y));
        }
    }
    else
    {
        for(i = 0; i < x->f_matrix_sizes.x; i++)
        {
            color_ref.lightness = x->f_lightness;
            for (j = 0; j < x->f_matrix_sizes.y; j++)
            {
                x->f_matrix_colorpanel[i][j] = color_ref;
                color_ref.lightness -= (1. / (float)(x->f_matrix_sizes.y - 1));
                if(x->f_matrix_colorpanel[i][j].lightness > 1)
                    x->f_matrix_colorpanel[i][j].lightness = 1. - (x->f_matrix_colorpanel[i][j].lightness - 1.);
                else if(x->f_matrix_colorpanel[i][j].lightness < 0.)
                    x->f_matrix_colorpanel[i][j].lightness = -x->f_matrix_colorpanel[i][j].lightness;
            }
            color_ref.hue += (1. / (float)(x->f_matrix_sizes.x));
        }
    }

    ebox_invalidate_layer((t_ebox *)x, cream_sym_background_layer);
    ebox_invalidate_layer((t_ebox *)x, cream_sym_hover_layer);
    ebox_invalidate_layer((t_ebox *)x, cream_sym_picked_layer);
    ebox_redraw((t_ebox *)x);
}

static t_pd_err colorpanel_matrix_set(t_colorpanel *x, t_object *attr, int ac, t_atom *av)
{
    int i;
    if(ac > 1 && av && atom_gettype(av) == A_FLOAT && atom_gettype(av+1) == A_FLOAT)
    {
        for(i = 0; i < x->f_matrix_sizes.x; i++)
        {
            free(x->f_matrix_colorpanel[i]);
        }
        free(x->f_matrix_colorpanel);

        x->f_matrix_sizes.x = (int)pd_clip_min(atom_getfloat(av), 1);
        x->f_matrix_sizes.y = (int)pd_clip_min(atom_getfloat(av+1), 1);
        x->f_matrix_colorpanel   = (t_hsla **)malloc(x->f_matrix_sizes.x * sizeof(t_hsla *));
        for(i = 0; i < x->f_matrix_sizes.x; i++)
        {
            x->f_matrix_colorpanel[i]   = (t_hsla *)malloc(x->f_matrix_sizes.y * sizeof(t_hsla));
        }
        colorpanel_computecolors(x);
    }
    return 0;
}

static t_pd_err colorpanel_saturation_set(t_colorpanel *x, t_object *attr, int ac, t_atom *av)
{
    if(ac && av && atom_gettype(av) == A_FLOAT)
    {
        x->f_saturation = pd_clip(atom_getfloat(av), 0., 1.);
        colorpanel_computecolors(x);
    }
    return 0;
}

static t_pd_err colorpanel_hue_set(t_colorpanel *x, t_object *attr, int ac, t_atom *av)
{
    if(ac && av && atom_gettype(av) == A_FLOAT)
    {
        x->f_hue = pd_clip(atom_getfloat(av), 0., 1.);
        colorpanel_computecolors(x);
    }
    return 0;
}

static t_pd_err colorpanel_lightness_set(t_colorpanel *x, t_object *attr, int ac, t_atom *av)
{
    if(ac && av && atom_gettype(av) == A_FLOAT)
    {
        x->f_lightness = pd_clip(atom_getfloat(av), 0., 1.);
        colorpanel_computecolors(x);
    }
    return 0;
}

static t_pd_err colorpanel_reverse_set(t_colorpanel *x, t_object *attr, int ac, t_atom *av)
{
    if(ac && av && atom_gettype(av) == A_FLOAT)
    {
        x->f_reverse = pd_clip((int)atom_getfloat(av), 0., 1.);
        colorpanel_computecolors(x);
    }
    return 0;
}

extern "C" void setup_c0x2ecolorpanel(void)
{
    t_eclass *c;
    
    c = eclass_new("c.colorpanel", (method)colorpanel_new, (method)colorpanel_free, (short)sizeof(t_colorpanel), 0L, A_GIMME, 0);
    eclass_guiinit(c, 0);
    
    eclass_addmethod(c, (method) colorpanel_paint,           "paint",            A_NULL, 0);
    eclass_addmethod(c, (method) colorpanel_notify,          "notify",           A_NULL, 0);
    eclass_addmethod(c, (method) colorpanel_getdrawparams,   "getdrawparams",    A_NULL, 0);
    eclass_addmethod(c, (method) colorpanel_oksize,          "oksize",           A_NULL, 0);
    eclass_addmethod(c, (method) colorpanel_set,             "set",              A_GIMME,0);
    eclass_addmethod(c, (method) colorpanel_list,            "list",             A_GIMME,0);
    eclass_addmethod(c, (method) colorpanel_output,          "bang",             A_NULL, 0);
    
    eclass_addmethod(c, (method) colorpanel_mousemove,       "mousemove",        A_NULL, 0);
    eclass_addmethod(c, (method) colorpanel_mousedown,       "mousedown",        A_NULL, 0);
    eclass_addmethod(c, (method) colorpanel_mousedown,       "mousedrag",        A_NULL, 0);
    eclass_addmethod(c, (method) colorpanel_mouseleave,      "mouseleave",       A_NULL, 0);
    
    eclass_addmethod(c, (method) colorpanel_preset,          "preset",           A_NULL, 0);
    
    CLASS_ATTR_INVISIBLE            (c, "send", 1);
    CLASS_ATTR_DEFAULT              (c, "size", 0, "181 105");
    
    CLASS_ATTR_FLOAT_ARRAY          (c, "matrix", 0, t_colorpanel, f_matrix_sizes, 2);
    CLASS_ATTR_LABEL                (c, "matrix", 0, "Matrix Size");
    CLASS_ATTR_ACCESSORS			(c, "matrix", NULL, colorpanel_matrix_set);
    CLASS_ATTR_ORDER                (c, "matrix", 0, "1");
    CLASS_ATTR_DEFAULT              (c, "matrix", 0, "24. 13.");
    CLASS_ATTR_SAVE                 (c, "matrix", 0);
    
    CLASS_ATTR_CHAR                 (c, "reverse", 0, t_colorpanel, f_reverse);
    CLASS_ATTR_LABEL                (c, "reverse", 0, "Matrix Reversed");
    CLASS_ATTR_ACCESSORS			(c, "reverse", NULL, colorpanel_reverse_set);
    CLASS_ATTR_ORDER                (c, "reverse", 0, "1");
    CLASS_ATTR_DEFAULT              (c, "reverse", 0, "0");
    CLASS_ATTR_SAVE                 (c, "reverse", 0);
    CLASS_ATTR_STYLE                (c, "reverse", 0, "onoff");
    
    CLASS_ATTR_FLOAT                (c, "saturation", 0, t_colorpanel, f_saturation);
    CLASS_ATTR_LABEL                (c, "saturation", 0, "Saturation");
    CLASS_ATTR_ACCESSORS			(c, "saturation", NULL, colorpanel_saturation_set);
    CLASS_ATTR_ORDER                (c, "saturation", 0, "1");
    CLASS_ATTR_DEFAULT              (c, "saturation", 0, "1.");
    CLASS_ATTR_SAVE                 (c, "saturation", 0);
    CLASS_ATTR_STYLE                (c, "saturation", 0, "number");
    CLASS_ATTR_STEP                 (c, "saturation", 0.1);
    
    CLASS_ATTR_FLOAT                (c, "hue", 0, t_colorpanel, f_hue);
    CLASS_ATTR_LABEL                (c, "hue", 0, "Fist Hue");
    CLASS_ATTR_ACCESSORS			(c, "hue", NULL, colorpanel_hue_set);
    CLASS_ATTR_ORDER                (c, "hue", 0, "1");
    CLASS_ATTR_DEFAULT              (c, "hue", 0, "0.");
    CLASS_ATTR_SAVE                 (c, "hue", 0);
    CLASS_ATTR_STYLE                (c, "hue", 0, "number");
    CLASS_ATTR_STEP                 (c, "hue", 0.1);
    
    CLASS_ATTR_FLOAT                (c, "lightness", 0, t_colorpanel, f_lightness);
    CLASS_ATTR_LABEL                (c, "lightness", 0, "First Lightness");
    CLASS_ATTR_ACCESSORS			(c, "lightness", NULL, colorpanel_lightness_set);
    CLASS_ATTR_ORDER                (c, "lightness", 0, "1");
    CLASS_ATTR_DEFAULT              (c, "lightness", 0, "1.");
    CLASS_ATTR_SAVE                 (c, "lightness", 0);
    CLASS_ATTR_STYLE                (c, "lightness", 0, "number");
    CLASS_ATTR_STEP                 (c, "lightness", 0.1);
    
    CLASS_ATTR_RGBA                 (c, "bgcolor", 0, t_colorpanel, f_color_background);
    CLASS_ATTR_LABEL                (c, "bgcolor", 0, "Background Color");
    CLASS_ATTR_ORDER                (c, "bgcolor", 0, "1");
    CLASS_ATTR_DEFAULT_SAVE_PAINT   (c, "bgcolor", 0, "0.75 0.75 0.75 1.");
    CLASS_ATTR_STYLE                (c, "bgcolor", 0, "color");
    
    CLASS_ATTR_RGBA                 (c, "bdcolor", 0, t_colorpanel, f_color_border);
    CLASS_ATTR_LABEL                (c, "bdcolor", 0, "Border Color");
    CLASS_ATTR_ORDER                (c, "bdcolor", 0, "2");
    CLASS_ATTR_DEFAULT_SAVE_PAINT   (c, "bdcolor", 0, "0.5 0.5 0.5 1.");
    CLASS_ATTR_STYLE                (c, "bdcolor", 0, "color");
    
    eclass_register(CLASS_BOX, c);
    colorpanel_class = c;
}





