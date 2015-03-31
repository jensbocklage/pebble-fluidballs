/* fluidballs, Copyright (c) 2000 by Peter Birtles <peter@bqdesign.com.au>
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *
 * Ported to X11 and xscreensaver by jwz, 27-Feb-2002.
 *
 * http://astronomy.swin.edu.au/~pbourke/modelling/fluid/
 *
 * Some physics improvements by Steven Barker <steve@blckknght.org>
 *
 * Port for pebble by Marcus Fritzsch
 */

/* Future ideas:
 * Specifying a distribution in the ball sizes (with a gamma curve, possibly).
 * Brownian motion, for that extra touch of realism.
 *
 * It would be nice to detect when there are more balls than fit in
 * the window, and scale the number of balls back.  Useful for the
 * xscreensaver-demo preview, which is often too tight by default.
 */

#include <pebble.h>

#if 0
#undef APP_LOG
#define APP_LOG(...)
#define START_TIME_MEASURE()
#define END_TIME_MEASURE()
#else
static unsigned int get_time(void)
{
   time_t s;
   uint16_t ms;
   time_ms(&s, &ms);
   return (s & 0xfffff) * 1000 + ms;
}

#define START_TIME_MEASURE() unsigned tm_0 = get_time()
#define END_TIME_MEASURE()                                                  \
   do                                                                       \
   {                                                                        \
      unsigned tm_1 = get_time();                                           \
      APP_LOG(APP_LOG_LEVEL_DEBUG, "%s: took %dms", __func__, tm_1 - tm_0); \
   } while (0)
#endif

#define M_PI 3.14159265358979323846

static float sqrtf(float f)
{
   float v = f * 0.5;
#define IT() v = (v + f / v) * 0.5f
   IT();
   IT();
   IT();
   IT();
   IT();  // 5 is not enough for a nice animation
   IT();
   IT();
   IT();  // 8 looks good
   IT();
   IT();  // 10 looks even better
#undef IT
   return v;
}

extern inline unsigned int GameRand(void)
{
   static unsigned int low = 16180, high = 31415;
   high = (high << 16) + (high >> 16);
   high += low;
   low += high;
   return high;
}

extern inline float frand(float m)
{
   float x = GameRand() / (float)-1u;
   return x * m;
}

#define NUMBALLS 50
#define GRAV 0.2

typedef struct
{
   GRect bounds;
   Window *window;
   Animation *anim;
   int count;                        /* number of balls */
   float accx;                       /* horizontal acceleration (wind) */
   float accy;                       /* vertical acceleration (gravity) */
   float vx[NUMBALLS], vy[NUMBALLS]; /* current ball velocities */
   float px[NUMBALLS], py[NUMBALLS]; /* current ball positions */
   float r[NUMBALLS];                /* ball radiuses */
   float m[NUMBALLS];                /* ball mass, precalculated */
   float e;                          /* coeficient of elasticity */
   float max_radius;                 /* largest radius of any ball */
} b_state;

static b_state s_state;

static void fluidballs_init(void)
{
   s_state.max_radius = 10;
   s_state.accx = 0;
   s_state.accy = GRAV;
   s_state.e = 0.97;
   s_state.count = NUMBALLS;

   for (int i = 0, end = s_state.count; i < end; i++)
   {
      float r = s_state.r[i] =
         (frand(s_state.max_radius * 0.75) + s_state.max_radius * 0.25);
      s_state.px[i] = frand(s_state.bounds.size.w - 2 * r) + r;
      s_state.py[i] = frand(s_state.bounds.size.h - 2 * r) + r;
      s_state.vx[i] = 0;  // frand(5) - 2.5;
      s_state.vy[i] = 0;  // frand(5) - 2.5;
      s_state.m[i] = r * r * r * M_PI * 4 / 3;

      APP_LOG(APP_LOG_LEVEL_DEBUG,
              "created ball %d: p=(%d, %d), v=(%d, %d), r=%d, m=%d", i,
              (int)s_state.px[i], (int)s_state.py[i], (int)s_state.vx[i],
              (int)s_state.vy[i], (int)s_state.r[i], (int)s_state.m[i]);
   }
}

/* Implements the laws of physics: move balls to their new positions.
 */
static void update_balls(void)
{
   APP_LOG(APP_LOG_LEVEL_DEBUG, "update_balls");

   float e = s_state.e;

   uint16_t collision_count = 0;
   START_TIME_MEASURE();

   /* For each ball, compute the influence of every other ball. */
   for (int a = 0, end = s_state.count; a < end - 1; a++)
   {
      float pxa = s_state.px[a], pya = s_state.py[a], ra = s_state.r[a],
            ma = s_state.m[a], vxa = s_state.vx[a], vya = s_state.vy[a];

      for (int b = a + 1; b < end; b++)
      {
         float pxb = s_state.px[b], pyb = s_state.py[b], rb = s_state.r[b];
         float d = (pxa - pxb) * (pxa - pxb) + (pya - pyb) * (pya - pyb);
         float dee2 = (ra + rb) * (ra + rb);

         if (d < dee2)
         {
            float mb = s_state.m[b];

            float vxb = s_state.vx[b];
            float vyb = s_state.vy[b];

            collision_count++;
            d = sqrtf(d);
            float rd = 1.f / d;
            float dd = ra + rb - d;

            float cdx = (pxb - pxa) * rd;
            float cdy = (pyb - pya) * rd;

            /* Move each ball apart from the other by half the
             * 'collision' distance.
             */
            pxa -= 0.5 * dd * cdx;
            pya -= 0.5 * dd * cdy;
            s_state.px[b] += 0.5 * dd * cdx;
            s_state.py[b] += 0.5 * dd * cdy;

            float vca =
               vxa * cdx + vya * cdy; /* the component of each velocity */
            float vcb =
               vxb * cdx + vyb * cdy; /* along the axis of the collision */

            /* elastic collison */
            float dva = (vca * (ma - mb) + vcb * 2 * mb) / (ma + mb) - vca;
            float dvb = (vcb * (mb - ma) + vca * 2 * ma) / (ma + mb) - vcb;

            dva *= e; /* some energy lost to inelasticity */
            dvb *= e;

#if 0
            dva += (frand (50) - 25) / ma;   /* q: why are elves so chaotic? */
            dvb += (frand (50) - 25) / mb;   /* a: brownian motion. */
#endif

            vxa += dva * cdx;
            vya += dva * cdy;
            vxb += dvb * cdx;
            vyb += dvb * cdy;

            s_state.vx[b] = vxb;
            s_state.vy[b] = vyb;
         }

         s_state.px[a] = pxa;
         s_state.py[a] = pya;
         s_state.vx[a] = vxa;
         s_state.vy[a] = vya;
      }
   }

   /* Force all balls to be on screen.
    */
   for (int a = 0, end = s_state.count; a < end; a++)
   {
      float r = s_state.r[a];
      if (s_state.px[a] < r)
      {
         s_state.px[a] = r;
         s_state.vx[a] = -s_state.vx[a] * e;
      }
      if (s_state.px[a] + r > s_state.bounds.size.w)
      {
         s_state.px[a] = s_state.bounds.size.w - r;
         s_state.vx[a] = -s_state.vx[a] * e;
      }
      if (s_state.py[a] < r)
      {
         s_state.py[a] = r;
         s_state.vy[a] = -s_state.vy[a] * e;
      }
      if (s_state.py[a] + r > s_state.bounds.size.h)
      {
         s_state.py[a] = s_state.bounds.size.h - r;
         s_state.vy[a] = -s_state.vy[a] * e;
      }
   }

   /* Apply gravity to all balls.
    */
   for (int a = 0; a < s_state.count; a++)
   {
      s_state.vx[a] += s_state.accx;
      s_state.vy[a] += s_state.accy;
      s_state.px[a] += s_state.vx[a];
      s_state.py[a] += s_state.vy[a];
   }

   END_TIME_MEASURE();

   APP_LOG(APP_LOG_LEVEL_DEBUG, "update_balls, had %d collisions",
           collision_count);
}

/* Erases the balls at their previous positions, and draws the new ones.
 */
static void repaint_balls(Layer *layer, GContext *ctx)
{
   APP_LOG(APP_LOG_LEVEL_DEBUG, "repaint_balls");

   START_TIME_MEASURE();

   graphics_context_set_antialiased(ctx, false);

   graphics_context_set_fill_color(ctx, GColorBlack);
   graphics_context_set_stroke_width(ctx, 0);

   graphics_fill_rect(ctx, s_state.bounds, 0, 0);

   graphics_context_set_fill_color(ctx, GColorWhite);
   for (int a = 0; a < s_state.count; a++)
   {
      graphics_fill_circle(
         ctx, (GPoint){.x = s_state.px[a], .y = s_state.py[a]}, s_state.r[a]);
   }

   END_TIME_MEASURE();
}

static void window_load(Window *window)
{
   layer_set_update_proc(window_get_root_layer(window), repaint_balls);
   animation_schedule(s_state.anim);
   APP_LOG(APP_LOG_LEVEL_DEBUG, "window %p loaded", window);
}

static void update_gravity(void)
{
#if 0
// untested
   AccelData adata;
   int e;
   if ((e = accel_service_peek(&adata)) < 0) {
      APP_LOG(APP_LOG_LEVEL_DEBUG, "Could not get accel data: %d", e);
      return;
   }
   s_state.accx = adata.x * 0.1;
   s_state.accy = adata.y * 0.1;
#elif 1
   static int u = 0;

   u++;

   const int frames = 120;

   switch (u / frames)
   {
   case 0:
      s_state.accx = 0;
      s_state.accy = GRAV;
      break;
   case 1:
      s_state.accx = GRAV;
      s_state.accy = 0;
      break;
   case 2:
      s_state.accx = 0;
      s_state.accy = -GRAV;
      break;
   case 3:
      s_state.accx = -GRAV;
      s_state.accy = 0;
      break;
   case 4:
      s_state.accx = 0;
      break;
   }

   // 6, let no grav last 2x as long
   if (u >= frames * 6)
   {
      u = 0;
   }
#endif
}

static void window_unload(Window *window)
{
   layer_set_update_proc(window_get_root_layer(window), NULL);
   animation_unschedule_all();
   APP_LOG(APP_LOG_LEVEL_DEBUG, "window %p unloaded", window);
}

static void anim_setup(Animation *anim) {}

static void anim_update(Animation *anim, const uint32_t d)
{
   layer_mark_dirty(window_get_root_layer(s_state.window));
   update_gravity();
   update_balls();
}

static void anim_teardown(Animation *anim) {}

static AnimationImplementation anim_impl = {
   .setup = anim_setup, .update = anim_update, .teardown = anim_teardown};

static void init(void)
{
   s_state.window = window_create();
   window_set_fullscreen(s_state.window, true);
   s_state.bounds = (GRect){.size = (GSize){.w = 144, .h = 168}};
   s_state.anim = animation_create();
   fluidballs_init();
   animation_set_duration(s_state.anim, ANIMATION_DURATION_INFINITE);
   animation_set_implementation(s_state.anim, &anim_impl);
   window_set_window_handlers(s_state.window,
                              (WindowHandlers){
                                 .load = window_load, .unload = window_unload,
                              });
   window_stack_push(s_state.window, false);
}

static void deinit(void) { window_destroy(s_state.window); }

int main(void)
{
   init();

   APP_LOG(APP_LOG_LEVEL_DEBUG, "Done initializing, pushed window: %p",
           s_state.window);

   app_event_loop();
   deinit();
}
