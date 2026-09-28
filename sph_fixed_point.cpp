#include <cinttypes>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;


#define WIDTH               128
#define HEIGHT              64

#define DISPLAY_WIDTH       32
#define DISPLAY_HEIGHT      16

#define DT                  0.05
#define SUBSTEPS            4
#define H                   DT / SUBSTEPS

#define NUM_PARTICLES       512
#define PARTICLES_PER_CELL  16

#define KERNEL_RADIUS       4
#define WALL_RESTITUTION    0.85
#define PRESSURE_MULT       40
#define VISCOSITY_MULT      0
#define DAMPING             0.0

#define CONTAINER_HYP       sqrt(WIDTH * WIDTH + HEIGHT * HEIGHT)
#define SPATIAL_GRID_SIZE   (int)ceil(CONTAINER_HYP / KERNEL_RADIUS)


int16_t gravityX = 0;
int16_t gravityY = -256;
int16_t currentRotation = 0;
int16_t lastRotation = currentRotation;


struct SpatialGridCell
{
    int array[PARTICLES_PER_CELL];
    uint16_t size;
};


int16_t smoothingKernel(int16_t distance)
{
    int16_t normalized = distance / KERNEL_RADIUS;
    int16_t value = max(1.f - normalized * normalized, 0);
    return value * value * value;
}

int16_t d_smoothingKernel(int16_t distance)
{
    int16_t normalized = distance / KERNEL_RADIUS;
    int16_t value = max(1 - normalized * normalized, 0);
    return -6 * value * value * normalized / KERNEL_RADIUS;
}


void spatialGridCoord(
    int16_t *predPosXs,
    int16_t *predPosYs,
    int index,
    int16_t *out_xi,
    int16_t *out_yi)
{
    int16_t xi =
        (int16_t)floor(
            (predPosXs[index] + CONTAINER_HYP * 0.5) / KERNEL_RADIUS);
    int16_t yi =
        (int16_t)floor(
            (predPosYs[index] + CONTAINER_HYP * 0.5) / KERNEL_RADIUS);

    if (xi < 0) xi = 0;
    if (xi >= SPATIAL_GRID_SIZE) xi = SPATIAL_GRID_SIZE - 1;

    if (yi < 0) yi = 0;
    if (yi >= SPATIAL_GRID_SIZE) yi = SPATIAL_GRID_SIZE - 1;

    *out_xi = xi;
    *out_yi = yi;
}

void updateSpatialGrid(
    int16_t *predPosXs,
    int16_t *predPosYs,
    SpatialGridCell *spatialGrid)
{
    for (int i = 0; i < SPATIAL_GRID_SIZE; i++)
    {
        for (int j = 0; j < SPATIAL_GRID_SIZE; j++)
        {
            spatialGrid[j * SPATIAL_GRID_SIZE + i].size = 0;
        }
    }

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        int16_t xi, yi;

        spatialGridCoord(
            predPosXs,
            predPosYs,
            i,
            &xi,
            &yi
        );

        if (spatialGrid[yi * SPATIAL_GRID_SIZE + xi].size == PARTICLES_PER_CELL - 1) continue;

        spatialGrid[yi * SPATIAL_GRID_SIZE + xi].array[
            spatialGrid[yi * SPATIAL_GRID_SIZE + xi].size
        ] = i;

        spatialGrid[yi * SPATIAL_GRID_SIZE + xi].size += 1;
    }
};

void nearbyIndices(
    int16_t *predPosXs,
    int16_t *predPosYs,
    SpatialGridCell *spatialGrid,
    int index,
    int *out_indices,
    int *out_indicesSize)
{
    int16_t xi, yi;

    spatialGridCoord(
        predPosXs,
        predPosYs,
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
            SpatialGridCell cell = spatialGrid[yj * SPATIAL_GRID_SIZE + xj];

            for (int i = 0; i < cell.size; i++)
            {
                out_indices[*out_indicesSize] = cell.array[i];
                (*out_indicesSize)++;
            }
        }
    }
}


int16_t calculateDensity(
    int16_t *predPosXs,
    int16_t *predPosYs,
    SpatialGridCell *spatialGrid,
    int index)
{
    int16_t res = 0;

    int indices[9 * PARTICLES_PER_CELL];
    int indicesSize = 0;

    nearbyIndices(
        predPosXs,
        predPosYs,
        spatialGrid,
        index,
        indices,
        &indicesSize
    );

    for (int j = 0; j < indicesSize; j++)
    {
        int i = indices[j];
        if (index == i) continue;

        int16_t dx = predPosXs[i] - predPosXs[index];
        int16_t dy = predPosYs[i] - predPosYs[index];
        int16_t dist = sqrt(dx * dx + dy * dy);

        if (dist > KERNEL_RADIUS) continue;

        res += smoothingKernel(dist);
    }

    return res;
}


void calculateForce(
    int16_t *predPosXs,
    int16_t *predPosYs,
    int16_t *velXs,
    int16_t *velYs,
    int16_t *densities,
    SpatialGridCell *spatialGrid,
    int index,
    int16_t *out_forceX,
    int16_t *out_forceY)
{
    int16_t forceX = gravityX;
    int16_t forceY = gravityY;

    int indices[9 * PARTICLES_PER_CELL];
    int indicesSize = 0;

    nearbyIndices(
        predPosXs,
        predPosYs,
        spatialGrid,
        index,
        indices,
        &indicesSize
    );

    for (int j = 0; j < indicesSize; j++)
    {
        int i = indices[j];
        if (index == i) continue;

        int16_t dx = predPosXs[i] - predPosXs[index];
        int16_t dy = predPosYs[i] - predPosYs[index];
        int16_t dist = sqrt(dx * dx + dy * dy);

        if (dist > KERNEL_RADIUS) continue;

        int16_t dirx = dist > 0 ? dx / dist : 0;
        int16_t diry = dist > 0 ? dy / dist : 0;

        int16_t density = densities[i];
        int16_t pressureSlope = d_smoothingKernel(dist);
        int16_t pressure = density * pressureSlope;

        forceX += pressure * dirx * PRESSURE_MULT;
        forceY += pressure * diry * PRESSURE_MULT;

        // viscosity
        int16_t dvx = velXs[i] - velXs[index];
        int16_t dvy = velYs[i] - velYs[index];

        int16_t viscosityInfluence = smoothingKernel(dist);

        forceX += dvx * viscosityInfluence * VISCOSITY_MULT;
        forceY += dvy * viscosityInfluence * VISCOSITY_MULT;
    }

    *out_forceX = forceX;
    *out_forceY = forceY;
}

void updateSim(
    int16_t *posXs,
    int16_t *posYs,
    int16_t *predPosXs,
    int16_t *predPosYs,
    int16_t *velXs,
    int16_t *velYs,
    int16_t *densities,
    SpatialGridCell *spatialGrid)
{
    int16_t dRotation = currentRotation - lastRotation;
    lastRotation = currentRotation;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        predPosXs[i] = posXs[i] + velXs[i] * H;
        predPosYs[i] = posYs[i] + velYs[i] * H;
    }


    updateSpatialGrid(predPosXs, predPosYs, spatialGrid);


    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        densities[i] = calculateDensity(
            predPosXs,
            predPosYs,
            spatialGrid,
            i
        );
    }


    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        int16_t force_x, force_y;

        calculateForce(
            predPosXs,
            predPosYs,
            velXs,
            velYs,
            densities,
            spatialGrid,
            i,
            &force_x,
            &force_y
        );

        velXs[i] += force_x * H;
        velYs[i] += force_y * H;

        velXs[i] *= (1 - DAMPING * H);
        velYs[i] *= (1 - DAMPING * H);

        posXs[i] += velXs[i] * H;
        posYs[i] += velYs[i] * H;
    }


    int16_t c = cos(currentRotation);
    int16_t s = sin(currentRotation);

    int16_t cd = cos(dRotation);
    int16_t sd = sin(dRotation);

    int16_t normalXs[] = {c, s, -c, -s};
    int16_t normalYs[] = {s, -c, -s, c};
    int16_t distances[] = {0.5f * WIDTH, 0.5f * HEIGHT, 0.5f * WIDTH, 0.5f * HEIGHT};

    // return;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        for (int w = 0; w < 4; w++)
        {
            int16_t normalX = normalXs[w];
            int16_t normalY = normalYs[w];

            int16_t distance = distances[w];
            int16_t pdot =
                normalX * (posXs[i] + normalX * distance) +
                normalY * (posYs[i] + normalY * distance);

            if (pdot < 0)
            {
                int16_t wallVelX = ((cd - 1) * posXs[i] - sd * posYs[i]) / H;
                int16_t wallVelY = (sd * posXs[i] + (cd - 1) * posYs[i]) / H;

                int16_t relVelX = velXs[i] - wallVelX;
                int16_t relVelY = velYs[i] - wallVelY;

                int16_t vdot = normalX * relVelX + normalY * relVelY;

                int16_t reflectedVelX =
                    relVelX - (1 + WALL_RESTITUTION) * vdot * normalX + wallVelX;
                int16_t reflectedVelY =
                    relVelY - (1 + WALL_RESTITUTION) * vdot * normalY + wallVelY;

                int16_t dtFrac = pdot / vdot;
                posXs[i] -= velXs[i] * dtFrac;
                posYs[i] -= velYs[i] * dtFrac;

                posXs[i] += reflectedVelX * dtFrac;
                posYs[i] += reflectedVelY * dtFrac;

                velXs[i] = reflectedVelX;
                velYs[i] = reflectedVelY;
            }
        }
    }
};




int main()
{
    int16_t posXs[NUM_PARTICLES] = {0};
    int16_t posYs[NUM_PARTICLES] = {0};

    int16_t counter = 0;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        int16_t x = fmod(counter, 1.f);
        int16_t y = (counter - x) * 0.03f;

        posXs[i] = (x - 0.5f) * WIDTH;
        posYs[i] = (y - 0.5f) * HEIGHT;
        counter += 0.03;
    }

    int16_t predPosXs[NUM_PARTICLES] = {0};
    int16_t predPosYs[NUM_PARTICLES] = {0};

    int16_t velXs[NUM_PARTICLES] = {0};
    int16_t velYs[NUM_PARTICLES] = {0};

    int16_t densities[NUM_PARTICLES] = {0};
    SpatialGridCell spatialGrid[SPATIAL_GRID_SIZE * SPATIAL_GRID_SIZE];


    while (true)
    {
        updateSim(
            posXs,
            posYs,
            predPosXs,
            predPosYs,
            velXs,
            velYs,
            densities,
            spatialGrid
        );

        bool buf[DISPLAY_WIDTH][DISPLAY_HEIGHT] = {false};

        int16_t c = cos(-currentRotation);
        int16_t s = sin(-currentRotation);

        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            int16_t x = c * posXs[i] - s * posYs[i];
            int16_t y = s * posXs[i] + c * posYs[i];

            int xi = (int)((0.5f + x / WIDTH) * DISPLAY_WIDTH);
            int yi = (int)((0.5f + y / HEIGHT) * DISPLAY_HEIGHT);

            xi = min(xi, DISPLAY_WIDTH - 1);
            yi = min(yi, DISPLAY_HEIGHT - 1);

            buf[xi][yi] = true;
        }

        stringstream ss;

        ss << "\n\n\n\n\n\n";

        for (int yi = 0; yi < DISPLAY_HEIGHT; yi++)
        {
            for (int xi = 0; xi < DISPLAY_WIDTH; xi++)
            {
                ss << (buf[xi][DISPLAY_HEIGHT - yi - 1] ? "**" : "  ");
            }
            ss << "\n";
        }

        for (int i = 0; i < DISPLAY_WIDTH; i++)
        {
            ss << "--";
        }

        cout << ss.str();
    }

    return 0;
}
