#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <cairo.h>
#include "lib/sph.h"
#include "lib/kernels.h"
#include "lib/integrators.h"
#include "lib/simulation.h"
#include "lib/storage.h"
#include "lib/memory.h"
#include "lib/colormap.h"

#ifndef M_PI
#define M_PI 3.1415926535
#endif

#ifndef M_TAU
#define M_TAU 6.283185307
#endif

void randomize_particles(void)
{
	srand(time(NULL));
	for(int i = 0; i < pcount; i++)
	{
		struct particle *particle = &particles[i];
		particle->basic.m = 1.0 / pcount;
		particle->basic.r[0] = 0.5 * (((double)(rand())) / RAND_MAX);
		particle->basic.r[1] = ((double)(rand())) / RAND_MAX;
		particle->basic.r[2] = 0.0;
		particle->basic.v[0] = 0.1;
		for(int j = 1; j < 3; j++)
		{
			particle->basic.v[j] = 0.0;
		}
	}
}

void update_particles_volrad(void)
{
	for(int i = 0; i < pcount; i++)
	{
		struct particle *particle = &particles[i];
		particle->volume = particle->basic.m / particle->basic.d;
		particle->radius = cbrt(0.75 * particle->volume / M_PI);
		particle->basic.r[2] = 0.0;
	}
}

void correct_particles(void)
{
	for(int i = 0; i < pcount; i++)
	{
		struct particle *particle = &particles[i];
		if(particle->basic.r[0] < 0.0)
		{
			particle->basic.r[0] = 0.0;
			for(int j = 0; j < 3; j++)
			{
				particle->basic.v[j] = 0.0;
			}
		}
		else if(particle->basic.r[0] > clength * ccount[0])
		{
			particle->basic.r[0] = clength * ccount[0];
			for(int j = 0; j < 3; j++)
			{
				particle->basic.v[j] = 0.0;
			}
		}
		if(particle->basic.r[1] < 0.0)
		{
			particle->basic.r[1] = 0.0;
			if(particle->basic.v[1] < 0.0)
			{
				particle->basic.v[1] = 0.0;
			}
		}
	}
}

double *acceleration_gravity(double r[3], double b[3])
{
	double ag[3] = {0.0, -0.1, 0.0};
	for(int j = 0; j < 3; j++)
	{
		b[j] = ag[j];
	}
	return b;
}

double *acceleration_circular(double r[3], double b[3])
{
	double ag[3] = {0.0, 0.0, 0.0};
	ag[0] -= 0.2 * (r[1] - 0.5);
	ag[1] += 0.2 * (r[0] - 0.5);
	for(int j = 0; j < 3; j++)
	{
		b[j] = ag[j];
	}
	return b;
}

int width = 1000;
int height = 1000;
double zoom = 0.8;

int main(void)
{
	printf("main program started\n");

	h = 0.04;

	pcount = 1000;
	allocate_particles();

	clength = 0.04;
	ccount[0] = 25;
	ccount[1] = 35;
	ccount[2] = 1;
	allocate_cells();

	n = 1;
	allocate_sparticles();

	k = 0.1;
	d0 = 10.0;
	u = 0.1;
	t = 0.0;
	dt = 0.01;
	nt = 0;
	acceleration_external = NULL;
	integrator = integrate_euler_explicit;
	initial_conditions = randomize_particles;
	boundary_conditions = correct_particles;
	acceleration_external = acceleration_gravity;

	simulate_particles(false);
	nt = 10;

	compute_densities();

	load_colormap();

	for(int f = 0; f < 100; f++)
	{
		update_particles_volrad();

		printf("frame:%d\n", f);

		cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
		cairo_t *context = cairo_create(surface);

		cairo_save(context);
		cairo_translate(context, 0.0, 0.5 * ((double) height));
		cairo_scale(context, 1.0, -1.0);
		cairo_translate(context, 0.0, -0.5 * ((double) height));

		cairo_rectangle(context, 0.0, 0.0, (double) width, (double) height);
		cairo_set_source_rgb(context, 0.0, 0.0, 0.0);
		cairo_fill(context);

		cairo_save(context);
		cairo_translate(context, 0.5 * ((double) width), 0.5 * ((double) height));
		cairo_scale(context, (double) width, (double) height);
		cairo_scale(context, zoom, zoom);
		cairo_translate(context, -0.5, -0.5);

		for(int i = 0; i < pcount; i++)
		{
			struct particle *particle = &particles[i];
			cairo_arc(context, particle->basic.r[0], particle->basic.r[1], 0.01, 0, M_TAU);
			double speed = 0.0;
			for(int j = 0; j < 3; j++)
			{
				speed += pow(particle->basic.v[j], 2);
			}
			speed = sqrt(speed);
			double value = 10 * speed;
			int index = index_colormap(value);
			float color[3] = {
					colors[index][0],
					colors[index][1],
					colors[index][2]
			};
			cairo_set_source_rgb(context, color[0], color[1], color[2]);
			cairo_fill_preserve(context);
			cairo_set_line_width(context, 0.003);
			cairo_set_source_rgb(context, 1.0, 1.0, 1.0);
			cairo_stroke(context);
		}

		for(int x = 0; x < ccount[0]; x++)
		{
			for(int y = 0; y < ccount[1]; y++)
			{
				for(int z = 0; z < ccount[2]; z++)
				{
					struct cell *cell = &cells[x][y][z];
					cairo_rectangle(context, cell->center[0] - 0.5 * cell->length, cell->center[1] - 0.5 * cell->length, cell->length, cell->length);
					cairo_set_line_width(context, 0.003);
					cairo_set_source_rgb(context, 1.0, 1.0, 1.0);
					cairo_stroke(context);
				}
			}
		}

		cairo_restore(context);
		cairo_restore(context);

		cairo_destroy(context);
		char filename[25];
		sprintf(filename, "out/img/%05d.png", f);
		cairo_surface_write_to_png(surface, filename);
		cairo_surface_destroy(surface);

		simulate_particles(true);
	}

	deallocate_sparticles();
	deallocate_cells();
	deallocate_particles();

	return EXIT_SUCCESS;
}
