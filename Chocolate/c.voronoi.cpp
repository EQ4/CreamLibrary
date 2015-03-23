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

#include "Voronoi.hpp"
using namespace voro;

typedef Voronoi<double>::Point Point;
typedef struct _voronoi
{
	t_ebox      j_box;
    t_outlet*   f_out;
    Voronoi<double>*    f_voro;
	t_rgba		f_color_background;
	t_rgba		f_color_border;
	t_rgba		f_color_voronoi;
    t_rect      f_rect;
    //t_pt        f_pt;
    long       f_index;
    void*       f_attr;
} t_voronoi;

t_eclass *voronoi_class;

void *voronoi_new(t_symbol *s, int argc, t_atom *argv);
void voronoi_free(t_voronoi *x);
void voronoi_assist(t_voronoi *x, void *b, long m, long a, char *s);

void voronoi_output(t_voronoi *x, t_symbol* s, long argc, t_atom* argv);
void voronoi_float(t_voronoi *x, float f);
void voronoi_symbol(t_voronoi *x, t_symbol* s);

t_pd_err voronoi_notify(t_voronoi *x, t_symbol *s, t_symbol *msg, void *sender, void *data);

void voronoi_getdrawparams(t_voronoi *x, t_object *patcherview, t_edrawparams *params);
void voronoi_oksize(t_voronoi *x, t_rect *newrect);

void voronoi_paint(t_voronoi *x, t_object *view);
void draw_background(t_voronoi *x,  t_object *view, t_rect *rect);

void voronoi_mousedown(t_voronoi *x, t_object *patcherview, t_pt pt, long modifiers);
void voronoi_mouseup(t_voronoi *x, t_object *patcherview, t_pt pt, long modifiers);
void voronoi_mousemove(t_voronoi *x, t_object *patcherview, t_pt pt, long modifiers);

void voronoi_setter(t_object *x, t_eattr* attr, long argc, t_atom* argv)
{
    post(attr->name->s_name);
    post(attr->type->s_name);
}

void voronoi_getter(t_object *x, t_eattr* attr, long* argc, t_atom** argv)
{
    post(attr->name->s_name);
    post(attr->type->s_name);
    argc[0] = 1;
    argv[0] = (t_atom *)malloc(sizeof(t_atom));
    atom_setlong(*argv, 7);
}

extern "C" void setup_c0x2evoronoi(void)
{
	t_eclass *c;
    
	c = eclass_new("c.voronoi", (method)voronoi_new, (method)voronoi_free, (short)sizeof(t_voronoi), 0L, A_GIMME, 0);
    
	eclass_init(c, 0);
    cream_initclass(c);
    
	eclass_addmethod(c, (method) voronoi_assist,          "assist",           A_NULL, 0);
	eclass_addmethod(c, (method) voronoi_paint,           "paint",            A_NULL, 0);
	eclass_addmethod(c, (method) voronoi_notify,          "notify",           A_NULL, 0);
    eclass_addmethod(c, (method) voronoi_getdrawparams,   "getdrawparams",    A_NULL, 0);
    eclass_addmethod(c, (method) voronoi_oksize,          "oksize",           A_NULL, 0);
    eclass_addmethod(c, (method) voronoi_float,          "float",            A_FLOAT,0);
    eclass_addmethod(c, (method) voronoi_symbol,          "symbol",            A_SYMBOL,0);
    eclass_addmethod(c, (method) voronoi_output,          "bang",             A_NULL, 0);
    eclass_addmethod(c, (method) voronoi_output,          "list",             A_GIMME,0);
    //eclass_addmethod(c, (method) voronoi_output,          "anything",         A_GIMME,0);
    
    eclass_addmethod(c, (method) voronoi_mousedown,       "mousedown",        A_NULL, 0);
    eclass_addmethod(c, (method) voronoi_mousedown,       "mousedrag",        A_NULL, 0);
    eclass_addmethod(c, (method) voronoi_mouseup,         "mouseup",          A_NULL, 0);
    eclass_addmethod(c, (method) voronoi_mousemove,       "mousemove",          A_NULL, 0);
    
    CLASS_ATTR_INVISIBLE            (c, "fontname", 1);
    CLASS_ATTR_INVISIBLE            (c, "fontweight", 1);
    CLASS_ATTR_INVISIBLE            (c, "fontslant", 1);
    CLASS_ATTR_INVISIBLE            (c, "fontsize", 1);
	CLASS_ATTR_DEFAULT              (c, "size", 0, "16. 16.");
    
    CLASS_ATTR_LONG                 (c, "void", 0, t_voronoi, f_attr);
    CLASS_ATTR_ACCESSORS            (c, "void", voronoi_getter, voronoi_setter);
    
	CLASS_ATTR_RGBA                 (c, "bgcolor", 0, t_voronoi, f_color_background);
	CLASS_ATTR_LABEL                (c, "bgcolor", 0, "Background Color");
	CLASS_ATTR_ORDER                (c, "bgcolor", 0, "1");
	CLASS_ATTR_DEFAULT_SAVE_PAINT   (c, "bgcolor", 0, "0.75 0.75 0.75 1.");
	CLASS_ATTR_STYLE                (c, "bgcolor", 0, "color");
    
	CLASS_ATTR_RGBA                 (c, "bdcolor", 0, t_voronoi, f_color_border);
	CLASS_ATTR_LABEL                (c, "bdcolor", 0, "Border Color");
	CLASS_ATTR_ORDER                (c, "bdcolor", 0, "2");
	CLASS_ATTR_DEFAULT_SAVE_PAINT   (c, "bdcolor", 0, "0.5 0.5 0.5 1.");
	CLASS_ATTR_STYLE                (c, "bdcolor", 0, "color");
    
	CLASS_ATTR_RGBA                 (c, "bacolor", 0, t_voronoi, f_color_voronoi);
	CLASS_ATTR_LABEL                (c, "bacolor", 0, "Bang Color");
	CLASS_ATTR_ORDER                (c, "bacolor", 0, "3");
	CLASS_ATTR_DEFAULT_SAVE_PAINT   (c, "bacolor", 0, "0. 0. 0. 1.");
	CLASS_ATTR_STYLE                (c, "bacolor", 0, "color");
	
    eclass_register(CLASS_BOX, c);
	voronoi_class = c;
}

void *voronoi_new(t_symbol *s, int argc, t_atom *argv)
{
	t_voronoi *x =  NULL;
	t_binbuf* d;
    long flags;
	if (!(d = binbuf_via_atoms(argc,argv)))
		return NULL;
    
	x = (t_voronoi *)eobj_new(voronoi_class);
    flags = 0
    | EBOX_GROWLINK
    ;
	ebox_new((t_ebox *)x, flags);
    x->f_out = (t_outlet *)bangout((t_object *)x);
    x->f_voro = new Voronoi<double>();
    x->f_index = -1l;
	ebox_attrprocess_viabinbuf(x, d);
	ebox_ready((t_ebox *)x);

	return (x);
}

void voronoi_getdrawparams(t_voronoi *x, t_object *patcherview, t_edrawparams *params)
{
	params->d_borderthickness   = 2;
	params->d_cornersize        = 2;
    params->d_bordercolor       = x->f_color_border;
    params->d_boxfillcolor      = x->f_color_background;
}

void voronoi_oksize(t_voronoi *x, t_rect *newrect)
{
    newrect->width = pd_clip_min(newrect->width, 16.);
    newrect->height = pd_clip_min(newrect->height, 16.);
    if((int)newrect->width % 2 == 0)
        newrect->width++;
    if((int)newrect->height % 2 == 0)
        newrect->height++;
}

double myrand()
{
    return (double(rand()) / RAND_MAX) * 2. -1;
}

void voronoi_list(t_voronoi *x, t_symbol* s, long argc, t_atom* argv)
{
    if(argc && argv)
    {
        for(int i = 0; i < 100; i++)
        {
            Point p(myrand(), myrand(), myrand());
            p.normalize();
            x->f_voro->add(p);
        }

        x->f_voro->compute();
        ebox_invalidate_layer((t_ebox *)x, gensym("background_layer"));
        ebox_redraw((t_ebox *)x);
    }
    
}

void voronoi_float(t_voronoi *x, float f)
{
    
    x->f_voro->clear();
    if((int)f%2==0)
    {
        f = f+1;
    }
    for(int i = 0; i < f; i++)
    {
        Point p = Point::fromPolar(1., (double)i / (double)f * HOA_2PI, (double)i / (double)(f -1.) * HOA_PI);
        x->f_voro->add(p);
    }

    x->f_voro->compute();
    ebox_invalidate_layer((t_ebox *)x, gensym("background_layer"));
    ebox_redraw((t_ebox *)x);
}

void voronoi_symbol(t_voronoi *x, t_symbol* s)
{
    x->f_voro->clear();
    if(s == gensym("square"))
    {
        x->f_voro->add(Point(1., 1., 1.));
        x->f_voro->add(Point(-1., 1., 1.));
        x->f_voro->add(Point(-1., -1., 1.));
        x->f_voro->add(Point(1., -1., 1.));
        
        x->f_voro->add(Point(1., 1., -1.));
        x->f_voro->add(Point(-1., 1., -1.));
        x->f_voro->add(Point(-1., -1., -1.));
        x->f_voro->add(Point(1., -1., -1.));
    }
    else if(s == gensym("tethrahedron"))
    {
        const double oh = -(sqrt(2. / 3.) / sqrt(3. / 8.) - 1.);
        const double hc = sqrt(1. - oh * oh);
        const double el = asin(oh / sqrt(hc*hc + oh*oh));
        x->f_voro->add(Point::fromPolar(1., 0., HOA_PI2));
        x->f_voro->add(Point::fromPolar(1., 0., el));
        x->f_voro->add(Point::fromPolar(1., HOA_2PI / 3., el));
        x->f_voro->add(Point::fromPolar(1., 2. * HOA_2PI / 3., el));
    }
    else if(s == gensym("octahedron"))
    {
        x->f_voro->add(Point::fromPolar(1., 0., HOA_PI2));
        x->f_voro->add(Point::fromPolar(1., 0., 0.));
        x->f_voro->add(Point::fromPolar(1., HOA_PI2, 0.));
        x->f_voro->add(Point::fromPolar(1., 2. * HOA_PI2, 0.));
        x->f_voro->add(Point::fromPolar(1., 3. * HOA_PI2, 0.));
        
    }
    else if(s == gensym("icosahedron"))
    {
        x->f_voro->add(Point::fromPolar(1., 0., HOA_PI2));
        for(int i = 1; i < 6; i++)
        {
            x->f_voro->add(Point::fromPolar(1., double(i - 1.) / 5. * HOA_2PI, atan(0.5)));
            x->f_voro->add(Point::fromPolar(1., double(i - 1.) / 5. * HOA_2PI - HOA_PI / 5., -atan(0.5)));
        }
        x->f_voro->add(Point::fromPolar(1., 0., -HOA_PI2));
    }
    else if(s == gensym("dodecahedron"))
    {
        const double phi = (sqrt(5.) - 1.) / 2.; // The golden ratio
        const double R = 1. / sqrt(3.);
        const double a = R;
        const double b = R / phi;
        const double c = R * phi;
        
        for(int i = -1; i < 2; i += 2)
        {
            for(int j = -1; j < 2; j += 2)
            {
                x->f_voro->add(Point(0., i * c * R, -j * b * R));
                x->f_voro->add(Point(i * c * R, j * b * R, 0.));
                x->f_voro->add(Point(i * b * R, 0., -j * c * R));
                for (int k = -1; k < 2; k += 2)
                {
                    x->f_voro->add(Point(i * a * R, j * a * R, -k * a * R));
                }
            }
        }
    }
    
    x->f_voro->compute();
    ebox_invalidate_layer((t_ebox *)x, gensym("background_layer"));
    ebox_redraw((t_ebox *)x);
}

void voronoi_output(t_voronoi *x, t_symbol* s, long argc, t_atom* argv)
{
    x->f_voro->clear();
    for(int i = 0; i < 100; i++)
    {
        x->f_voro->add(Point(myrand(), myrand(), myrand()));
    }
    
    x->f_voro->compute();
    ebox_invalidate_layer((t_ebox *)x, gensym("background_layer"));
    ebox_redraw((t_ebox *)x);
}

void voronoi_free(t_voronoi *x)
{
	ebox_free((t_ebox *)x);
    delete x->f_voro;
}

void voronoi_assist(t_voronoi *x, void *b, long m, long a, char *s)
{
	;
}

t_pd_err voronoi_notify(t_voronoi *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
	if (msg == gensym("attr_modified"))
	{
		if(s == gensym("bgcolor") || s == gensym("bdcolor") || s == gensym("bacolor"))
		{
			ebox_invalidate_layer((t_ebox *)x, gensym("background_layer"));
		}
        ebox_redraw((t_ebox *)x);
	}
	return 0;
}

static const double k = 0.55228474983079356430692996582365594804286956787109;

static void rotate(const double cosz, const double sinz, t_pt& p1) noexcept
{
    const double rx = p1.x * cosz - p1.y * sinz;
    p1.y = p1.x * sinz + p1.y * cosz;
    p1.x = rx;
}

static void create_small_arc(const double r, const double start, const double extend, t_pt const& ct, t_pt& p2, t_pt& p3, t_pt& p4)
{
    const double a = extend;
    const double cosz = cos(a * 0.5 + start);
    const double sinz = sin(a * 0.5 + start);
    t_pt p1;
    p4.x = r * cos(a * 0.5);
    p4.y = r * sin(a * 0.5);
    p1.x = p4.x;
    p1.y = -p4.y;
    p2.x = p1.x + k * tan(a * 0.5) * p4.y;
    p2.y = p1.y + k * tan(a * 0.5) * p4.x;
    p3.x = p2.x;
    p3.y = -p2.y;
    
    rotate(cosz, sinz, p2); rotate(cosz, sinz, p3); rotate(cosz, sinz, p4);
    p2.x += ct.x; p2.y += ct.y; p3.x += ct.x; p3.y += ct.y; p4.x += ct.x; p4.y += ct.y;
}

static void graphics_arc_to(t_elayer *g, t_pt prev, float cx, float cy, float extend)
{
    t_pt p2, p3, p4, c = {cx, cy};
    double radius   = pd_radius(prev.x - cx, prev.y - cy);
    double angle    = pd_angle(prev.x - cx, prev.y - cy);
    
    while(extend > EPD_2PI)
    {
        extend -= EPD_2PI;
    }
    while(extend < -EPD_2PI)
    {
        extend += EPD_2PI;
    }
   
    while(fabs(extend) >= EPD_PI4)
    {
        if(extend < 0.)
        {
            create_small_arc(radius, angle, -EPD_PI4, c, p2, p3, p4);
            extend += EPD_PI4;
            angle  -= EPD_PI4;
        }
        else
        {
            create_small_arc(radius, angle, EPD_PI4, c, p2, p3, p4);
            extend -= EPD_PI4;
            angle  += EPD_PI4;
        }
        egraphics_curve_to(g, p2.x, p2.y, p3.x, p3.y,  p4.x, p4.y);
    }
    if(fabs(extend) > HOA_EPSILON)
    {
        create_small_arc(radius, angle, extend, c, p2, p3, p4);
        egraphics_curve_to(g, p2.x, p2.y, p3.x, p3.y,  p4.x, p4.y);
    }
}

void draw_arc(t_voronoi *x, t_object *view, t_rect *rect)
{
    t_matrix transform;
    t_elayer *g = ebox_start_layer((t_ebox *)x, gensym("arc_layer"), rect->width, rect->height);
    const double width = rect->width / 2.;
    if (g)
    {
        t_rgba black = {0., 0., 0., 1.};
        
        egraphics_matrix_init(&transform, 1, 0, 0, -1, width, width);
        egraphics_set_matrix(g, &transform);
        egraphics_set_color_rgba(g, &black);
        
     
        t_pt start = {0., width * 0.5};
        egraphics_move_to(g, start.x, start.y);
        egraphics_arc_to(g, 0., 0., HOA_PI);
        egraphics_stroke(g);

        ebox_end_layer((t_ebox*)x, gensym("arc_layer"));
    }
    ebox_paint_layer((t_ebox *)x, gensym("arc_layer"), 0., 0.);
}

void voronoi_paint(t_voronoi *x, t_object *view)
{
	ebox_get_rect_for_view((t_ebox *)x, &x->f_rect);
    draw_background(x, view, &x->f_rect);
    //draw_arc(x, view, &x->f_rect);
}

void draw_background(t_voronoi *x, t_object *view, t_rect *rect)
{
    t_matrix transform;
	t_elayer *g = ebox_start_layer((t_ebox *)x, gensym("background_layer"), rect->width, rect->height);
	if (g)
	{
        t_rgba red = {1., 0., 0., 1.};
        t_rgba blue = {0., 0., 1., 1.};
        t_rgba green = {0., 1., 0., 1.};
        const double width = rect->width / 2.;
        
        egraphics_matrix_init(&transform, 1, 0, 0, -1, rect->width / 2., rect->width / 2.);
        egraphics_set_matrix(g, &transform);
        
        egraphics_set_color_rgba(g, &x->f_color_border);
        egraphics_set_line_width(g, 1.);
        egraphics_circle(g, 0., 0., rect->width / 2.);
        egraphics_stroke(g);

        vector<Point>& t = x->f_voro->getPoints();
        
        for(int i = 0; i < t.size(); i++)
        {
            if(t[i].bounds.size() > 2)
            {
                float angle1 = pd_angle(t[i].bounds[0].x, t[i].bounds[0].y);
                float radius1= pd_radius(t[i].bounds[0].x, t[i].bounds[0].y);
                egraphics_move_to(g, t[i].bounds[0].x * width, t[i].bounds[0].y * width);
                for(int j = 1; j < t[i].bounds.size(); j++)
                {
                    const float angle2 = pd_angle(t[i].bounds[j].x, t[i].bounds[j].y);
                    const float radius2= pd_radius(t[i].bounds[j].x, t[i].bounds[j].y);
                    float extend = angle2 - angle1;
                    if(extend > HOA_PI)
                    {
                        extend -= HOA_2PI;
                    }
                    else if(extend < -HOA_PI)
                    {
                        extend += HOA_2PI;
                    }
                    
                    if(fabs(extend) > HOA_EPSILON  && fabs(radius2 - radius1) < HOA_EPSILON && fabs(fabs(extend) - HOA_PI)> 1e-6)
                    {
                        egraphics_arc_to(g, 0., 0., extend);
                    }
                    else
                    {
                        egraphics_line_to(g, t[i].bounds[j].x * width, t[i].bounds[j].y * width);
                    }
                    angle1 = angle2;
                    radius1 = radius2;
                }
                
                const float angle2 = pd_angle(t[i].bounds[0].x, t[i].bounds[0].y);
                const float radius2 = pd_radius(t[i].bounds[0].x, t[i].bounds[0].y);
                float extend = angle2 - angle1;
                if(extend > HOA_PI)
                {
                    extend -= HOA_2PI;
                }
                else if(extend < -HOA_PI)
                {
                    extend += HOA_2PI;
                }
           
                if(fabs(extend) > HOA_EPSILON  && fabs(radius2 - radius1) < HOA_EPSILON && fabs(fabs(extend) - HOA_PI)> 1e-6)
                {
                    egraphics_arc_to(g, 0., 0., extend);
                }
                else
                {
                    egraphics_line_to(g, t[i].bounds[0].x * width, t[i].bounds[0].y * width);
                }
                
                egraphics_set_color_rgba(g, &green);
                egraphics_fill_preserve(g);
                egraphics_set_color_rgba(g, &red);
                egraphics_stroke(g);
                /*
                egraphics_set_color_rgba(g, &blue);
                egraphics_circle(g, t[i].bounds[0].x * width, t[i].bounds[0].y * width, 3.);
                egraphics_fill(g);
                post("%f", t[i].bounds[0].z);
                for(int j = 1; j < t[i].bounds.size(); j++)
                {
                    egraphics_circle(g, t[i].bounds[j].x * width, t[i].bounds[j].y * width, 3.);
                    egraphics_fill(g);
                    post("%f", t[i].bounds[j].z);
                }*/
            }
        }
        
        egraphics_set_color_rgba(g, &red);
        for(int i = 0; i < t.size(); i++)
        {
            if(t[i].z >= 0 && t[i].bounds.size() > 2)
            {
                egraphics_circle(g, t[i].x  * rect->width / 2.,t[i].y  * rect->width / 2., 3.);
                egraphics_fill(g);
            }
            
        }
        ebox_end_layer((t_ebox*)x, gensym("background_layer"));
	}
	ebox_paint_layer((t_ebox *)x, gensym("background_layer"), 0., 0.);
}

void voronoi_mousedown(t_voronoi *x, t_object *patcherview, t_pt pt, long modifiers)
{
    /*
    const double w = x->f_rect.width * 0.5;
    Point p(pt.x / w - 1., -(pt.y / w - 1.), 0.);
    p = Point::fromPolar(1., p.azimuth(), (1. - p.radius()) * HOA_PI2);
    
    vector<Triangle> t = x->f_voro->getTriangles();
    if(t[0].d >= p.greatCircleDistance(t[0].p))
    {
        post("inside");
    }
    else
    {
        post("outside");
    }
    x->f_pt.x = p.x;
    x->f_pt.y = p.y;
    ebox_invalidate_layer((t_ebox *)x, gensym("background_layer"));
    ebox_redraw((t_ebox *)x);
     */
}

bool near(Point const& p1, Point const& p2)
{
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)) < 0.1;
}

void voronoi_mousemove(t_voronoi *x, t_object *patcherview, t_pt pt, long modifiers)
{
    const double w = x->f_rect.width * 0.5;
    Point p(pt.x / w - 1., -(pt.y / w - 1.), 0.);
    ulong index = x->f_index;
    x->f_index = -1l;
    vector<Point> const& t = x->f_voro->getPoints();
    for(int i = 0; i < t.size(); i++)
    {
        if(near(p, t[i]) && t[i].z >= 0)
        {
            x->f_index = i;
        }
    }
    
    if(x->f_index != index)
    {
        ebox_invalidate_layer((t_ebox *)x, gensym("background_layer"));
        ebox_redraw((t_ebox *)x);
    }
    
}

void voronoi_mouseup(t_voronoi *x, t_object *patcherview, t_pt pt, long modifiers)
{
    ;
}






