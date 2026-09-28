#include <cinttypes>
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>

using namespace std;


#define WIDTH               128
#define HEIGHT              64

#define DT                  0.008
#define SUBSTEPS            4
#define H                   DT / SUBSTEPS

#define NUM_PARTICLES       256
#define PARTICLES_PER_CELL  8

#define KERNEL_RADIUS       1
#define WALL_RESTITUTION    0.85
#define PRESSURE_MULT       20
#define VISCOSITY_MULT      0
#define DAMPING             0.0

#define CONTAILER_HYP       sqrt(WIDTH * WIDTH + HEIGHT * HEIGHT)
#define SPATIAL_GRID_SIZE   (int)ceil(CONTAILER_HYP / KERNEL_RADIUS)


random_device dev;
mt19937 rng(dev());
uniform_real_distribution<> uniform(0.f, 1.f);


float gravity_x = 0;
float gravity_y = -10;
float current_rotation = 0;
float last_rotation = current_rotation;
float rotation_damp = 0.3;


struct SpatialGridCell
{
    int array[PARTICLES_PER_CELL];
    uint16_t size;
};


float smoothing_kernel(float distance)
{
    float normalized = distance / KERNEL_RADIUS;
    float value = max(1.f - normalized * normalized, 0.f);
    return value * value * value;
}

float d_smoothing_kernel(float distance)
{
    float normalized = distance / KERNEL_RADIUS;
    float value = max(1 - normalized * normalized, 0.f);
    return -6.f * value * value * normalized / KERNEL_RADIUS;
}


void spatial_grid_coord(
    float *pred_pos_xs,
    float *pred_pos_ys,
    int index,
    int16_t *out_xi,
    int16_t *out_yi)
{
    int16_t xi =
        (int16_t)floor(
            (pred_pos_xs[index] + CONTAILER_HYP * 0.5) / KERNEL_RADIUS);
    int16_t yi =
        (int16_t)floor(
            (pred_pos_ys[index] + CONTAILER_HYP * 0.5) / KERNEL_RADIUS);

    if (xi < 0) xi = 0;
    if (xi >= SPATIAL_GRID_SIZE) xi = SPATIAL_GRID_SIZE - 1;

    if (yi < 0) yi = 0;
    if (yi >= SPATIAL_GRID_SIZE) yi = SPATIAL_GRID_SIZE - 1;

    *out_xi = xi;
    *out_yi = yi;
}

void update_spatial_grid(
    float *pred_pos_xs,
    float *pred_pos_ys,
    SpatialGridCell *spatial_grid)
{
    for (int i = 0; i < SPATIAL_GRID_SIZE; i++)
    {
        for (int j = 0; j < SPATIAL_GRID_SIZE; j++)
        {
            spatial_grid[j * SPATIAL_GRID_SIZE + i].size = 0;
        }
    }

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        int16_t xi, yi;

        spatial_grid_coord(
            pred_pos_xs,
            pred_pos_ys,
            i,
            &xi,
            &yi
        );

        spatial_grid[yi * SPATIAL_GRID_SIZE + xi].array[
            spatial_grid[yi * SPATIAL_GRID_SIZE + xi].size
        ] = i;

        spatial_grid[yi * SPATIAL_GRID_SIZE + xi].size += 1;
    }
};

void nearby_indices(
    float *pred_pos_xs,
    float *pred_pos_ys,
    SpatialGridCell *spatial_grid,
    int index,
    int *out_indices,
    int *out_indices_size)
{
    int16_t xi, yi;

    spatial_grid_coord(
        pred_pos_xs,
        pred_pos_ys,
        index,
        &xi,
        &yi
    );

    int xiMin = max(xi - 1, 0);
    int xiMax = min(xi + 1, SPATIAL_GRID_SIZE - 1);
    int yiMin = max(yi - 1, 0);
    int yiMax = min(yi + 1, SPATIAL_GRID_SIZE - 1);

    for (int xj = xiMin; xj <= xiMax; xj++)
    {
        for (int yj = yiMin; yj <= yiMax; yj++)
        {
            SpatialGridCell cell = spatial_grid[yj * SPATIAL_GRID_SIZE + xj];

            for (int i = 0; i < cell.size; i++)
            {
                out_indices[*out_indices_size] = cell.array[i];
                *out_indices_size++;
            }
        }
    }
}


float calculate_density(
    float *pred_pos_xs,
    float *pred_pos_ys,
    SpatialGridCell *spatial_grid,
    int index)
{
    float res = 0;

    // xi, yi = spatialGridCoord(index);
    int indices[9 * PARTICLES_PER_CELL];
    int indices_size = 0;

    nearby_indices(
        pred_pos_xs,
        pred_pos_ys,
        spatial_grid,
        index,
        indices,
        &indices_size
    );

    for (int j = 0; j < indices_size; j++)
        // for (int i = 0; i < NUM_PARTICLES; i++)
    {
        int i = indices[j];
        if (index == i) continue;

        float dx = pred_pos_xs[i] - pred_pos_xs[index];
        float dy = pred_pos_ys[i] - pred_pos_ys[index];
        float dist = sqrt(dx * dx + dy * dy);

        if (dist > KERNEL_RADIUS) continue;

        res += smoothing_kernel(dist);
    }

    return res;
}


void calculate_force(
    float *pred_pos_xs,
    float *pred_pos_ys,
    float *vel_xs,
    float *vel_ys,
    float *densities,
    int index,
    float *out_force_x,
    float *out_force_y)
{
    float resx = gravity_x;
    float resy = gravity_y;

    // xi, yi = spatialGridCoord(index);

    // for (i of nearby_indices(index))
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        if (index == i) continue;

        float dx = pred_pos_xs[i] - pred_pos_xs[index];
        float dy = pred_pos_ys[i] - pred_pos_ys[index];
        float dist = sqrt(dx * dx + dy * dy);

        if (dist > KERNEL_RADIUS) continue;

        // pressure
        // float dirx = dist > 0 ? dx / dist : random() - 0.5;
        // float diry = dist > 0 ? dy / dist : random() - 0.5;
        float dirx = dist > 0 ? dx / dist : 0;
        float diry = dist > 0 ? dy / dist : 0;

        float density = densities[i];
        float pressureSlope = d_smoothing_kernel(dist);
        float pressure = density * pressureSlope;

        resx += pressure * dirx * PRESSURE_MULT;
        resy += pressure * diry * PRESSURE_MULT;

        // viscosity
        float dvx = vel_xs[i] - vel_xs[index];
        float dvy = vel_ys[i] - vel_ys[index];

        float viscosityInfluence = smoothing_kernel(dist);

        resx += dvx * viscosityInfluence * VISCOSITY_MULT;
        resy += dvy * viscosityInfluence * VISCOSITY_MULT;
    }

    *out_force_x = resx;
    *out_force_y = resy;
}

void sim_update(
    float *pos_xs,
    float *pos_ys,
    float *pred_pos_xs,
    float *pred_pos_ys,
    float *vel_xs,
    float *vel_ys,
    float *densities,
    SpatialGridCell *spatial_grid)
{
    // current_rotation = last_rotation + (sliderRotation - last_rotation) * rotation_damp;
    // current_rotation = last_rotation + 0.8 * H;
    // current_rotation = -0.9

    float d_rotation = current_rotation - last_rotation;
    last_rotation = current_rotation;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        pred_pos_xs[i] = pos_xs[i] + vel_xs[i] * H;
        pred_pos_ys[i] = pos_ys[i] + vel_ys[i] * H;
    }


    update_spatial_grid(pred_pos_xs, pred_pos_ys, spatial_grid);


    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        densities[i] = calculate_density(
            pred_pos_xs,
            pred_pos_ys,
            spatial_grid,
            i
        );
    }



    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        float force_x, force_y;

        calculate_force(
            pred_pos_xs,
            pred_pos_ys,
            vel_xs,
            vel_ys,
            densities,
            i,
            &force_x,
            &force_y
        );

        vel_xs[i] += force_x * H;
        vel_ys[i] += force_y * H;

        vel_xs[i] *= (1 - DAMPING * H);
        vel_ys[i] *= (1 - DAMPING * H);

        pos_xs[i] += vel_xs[i] * H;
        pos_ys[i] += vel_ys[i] * H;
    }



    float c = cos(current_rotation);
    float s = sin(current_rotation);

    float cd = cos(d_rotation);
    float sd = sin(d_rotation);

    float normal_xs[] = {c, s, -c, -s};
    float normal_ys[] = {s, -c, -s, c};
    float distances[] = {0.5f * WIDTH, 0.5f * HEIGHT, 0.5f * WIDTH, 0.5f * HEIGHT};

    // return;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        for (int w = 0; w < 4; w++)
        {
            float normal_x = normal_xs[w];
            float normal_y = normal_ys[w];

            float distance = distances[w];
            float pdot =
                normal_x * (pos_xs[i] + normal_x * distance) +
                normal_y * (pos_ys[i] + normal_y * distance);

            if (pdot < 0)
            {
                float wall_vel_x = ((cd - 1) * pos_xs[i] - sd * pos_ys[i]) / H;
                float wall_vel_y = (sd * pos_xs[i] + (cd - 1) * pos_ys[i]) / H;

                float rel_vel_x = vel_xs[i] - wall_vel_x;
                float rel_vel_y = vel_ys[i] - wall_vel_y;

                float vdot = normal_x * rel_vel_x + normal_y * rel_vel_y;

                float reflected_vel_x =
                    rel_vel_x - (1 + WALL_RESTITUTION) * vdot * normal_x + wall_vel_x;
                float reflected_vel_y =
                    rel_vel_y - (1 + WALL_RESTITUTION) * vdot * normal_y + wall_vel_y;

                float dt_frac = pdot / vdot;
                pos_xs[i] -= vel_xs[i] * dt_frac;
                pos_ys[i] -= vel_ys[i] * dt_frac;

                pos_xs[i] += reflected_vel_x * dt_frac;
                pos_ys[i] += reflected_vel_y * dt_frac;

                vel_xs[i] = reflected_vel_x;
                vel_ys[i] = reflected_vel_y;
            }
        }
    }
};




int main()
{
    float pos_xs[NUM_PARTICLES] = {0};
    float pos_ys[NUM_PARTICLES] = {0};

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        pos_xs[i] = (uniform(rng) - 0.5f) * WIDTH;
        pos_ys[i] = (uniform(rng) - 0.5f) * HEIGHT;
    }

    float pred_pos_xs[NUM_PARTICLES] = {0};
    float pred_pos_ys[NUM_PARTICLES] = {0};

    float vel_xs[NUM_PARTICLES] = {0};
    float vel_ys[NUM_PARTICLES] = {0};

    float densities[NUM_PARTICLES] = {0};
    SpatialGridCell spatial_grid[SPATIAL_GRID_SIZE * SPATIAL_GRID_SIZE];


    const int DISPLAY_WIDTH = 32;
    const int DISPLAY_HEIGHT = 16;

    while (true)
    {
        sim_update(
            pos_xs,
            pos_ys,
            pred_pos_xs,
            pred_pos_ys,
            vel_xs,
            vel_ys,
            densities,
            spatial_grid
        );

        bool buf[DISPLAY_WIDTH + 1][DISPLAY_HEIGHT + 1] = {false};

        float c = cos(-current_rotation);
        float s = sin(-current_rotation);

        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            float x = c * pos_xs[i] - s * pos_ys[i];
            float y = s * pos_xs[i] + c * pos_ys[i];

            int xi = (int)((0.5f + x / WIDTH) * DISPLAY_WIDTH);
            int yi = (int)((0.5f + y / HEIGHT) * DISPLAY_HEIGHT);

            buf[xi][yi] = true;
        }

        stringstream ss;

        ss << "\n\n\n\n\n\n";

        for (int yi = 0; yi < DISPLAY_HEIGHT + 1; yi++)
        {
            for (int xi = 0; xi < DISPLAY_WIDTH + 1; xi++)
            {
                ss << (buf[xi][DISPLAY_WIDTH - yi] ? "XX" : "  ");
            }
            ss << "\n";
        }

        cout << ss.str();
    }

    return 0;
}
