#include <cinttypes>
#include <cmath>
#include <iostream>
#include <sstream>
#include "fixed16.cpp"

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
#define PRESSURE_MULT       20
#define VISCOSITY_MULT      0
#define DAMPING             0.0

#define CONTAINER_HYP       sqrt(WIDTH * WIDTH + HEIGHT * HEIGHT)
#define SPATIAL_GRID_SIZE   (int)ceil(CONTAINER_HYP / KERNEL_RADIUS)



const fixed16 h  = fixed16(static_cast<float>(H));

const fixed16 kernelRadius  = fixed16(static_cast<float>(KERNEL_RADIUS));
const fixed16 pressuerMult  = fixed16(static_cast<float>(PRESSURE_MULT));
const fixed16 viscosityMult = fixed16(static_cast<float>(VISCOSITY_MULT));


fixed16 gravityX = fixed16(0.f);
fixed16 gravityY = fixed16(-10.f);
fixed16 currentRotation = fixed16(0.f);
fixed16 lastRotation = currentRotation;


struct SpatialGridCell
{
    int array[PARTICLES_PER_CELL];
    uint16_t size;
};


fixed16 smoothingKernel(fixed16 distance)
{
    fixed16 normalized = distance / kernelRadius;
    fixed16 value = max(
        fixed16(1.f) - normalized * normalized,
        fixed16(0.f)
    );
    return value * value * value;
}

fixed16 d_smoothingKernel(fixed16 distance)
{
    fixed16 normalized = distance / kernelRadius;
    fixed16 value = max(
        fixed16(1.f) - normalized * normalized,
        fixed16(0.f)
    );
    return fixed16(-6.f) * value * value * normalized / kernelRadius;
}


void spatialGridCoord(
    fixed16 *predPosXs,
    fixed16 *predPosYs,
    int index,
    int16_t *out_xi,
    int16_t *out_yi)
{
    fixed16 offset = fixed16(static_cast<float>(CONTAINER_HYP * 0.5f));

    int16_t xi = floor((predPosXs[index] + offset) / kernelRadius).to_int();
    int16_t yi = floor((predPosYs[index] + offset) / kernelRadius).to_int();

    if (xi < 0) xi = 0;
    if (xi >= SPATIAL_GRID_SIZE) xi = SPATIAL_GRID_SIZE - 1;

    if (yi < 0) yi = 0;
    if (yi >= SPATIAL_GRID_SIZE) yi = SPATIAL_GRID_SIZE - 1;

    *out_xi = xi;
    *out_yi = yi;
}

void updateSpatialGrid(
    fixed16 *predPosXs,
    fixed16 *predPosYs,
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
    fixed16 *predPosXs,
    fixed16 *predPosYs,
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


fixed16 calculateDensity(
    fixed16 *predPosXs,
    fixed16 *predPosYs,
    SpatialGridCell *spatialGrid,
    int index)
{
    fixed16 res = fixed16(0.f);

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

        fixed16 dx = predPosXs[i] - predPosXs[index];
        fixed16 dy = predPosYs[i] - predPosYs[index];
        fixed16 dist = sqrt(dx * dx + dy * dy);

        if (dist > kernelRadius) continue;

        res += smoothingKernel(dist);
    }

    return res;
}


void calculateForce(
    fixed16 *predPosXs,
    fixed16 *predPosYs,
    fixed16 *velXs,
    fixed16 *velYs,
    fixed16 *densities,
    SpatialGridCell *spatialGrid,
    int index,
    fixed16 *out_forceX,
    fixed16 *out_forceY)
{
    fixed16 forceX = gravityX;
    fixed16 forceY = gravityY;

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

        fixed16 dx = predPosXs[i] - predPosXs[index];
        fixed16 dy = predPosYs[i] - predPosYs[index];
        fixed16 dist = sqrt(dx * dx + dy * dy);

        if (dist > kernelRadius) continue;

        fixed16 dirx = dist > fixed16(0.f)
            ? dx / dist
            : fixed16(0.f);
        fixed16 diry = dist > fixed16(0.f)
            ? dy / dist
            : fixed16(0.f);

        fixed16 density = densities[i];
        fixed16 pressureSlope = d_smoothingKernel(dist);
        fixed16 pressure = density * pressureSlope;

        forceX += pressure * dirx * pressuerMult;
        forceY+= pressure * diry * pressuerMult;

        // viscosity
        fixed16 dvx = velXs[i] - velXs[index];
        fixed16 dvy = velYs[i] - velYs[index];

        fixed16 viscosityInfluence = smoothingKernel(dist);

        forceX += dvx * viscosityInfluence * viscosityMult;
        forceY+= dvy * viscosityInfluence * viscosityMult;
    }

    *out_forceX = forceX;
    *out_forceY = forceY;
}

void updateSim(
    fixed16 *posXs,
    fixed16 *posYs,
    fixed16 *predPosXs,
    fixed16 *predPosYs,
    fixed16 *velXs,
    fixed16 *velYs,
    fixed16 *densities,
    SpatialGridCell *spatialGrid)
{
    fixed16 dRotation = currentRotation - lastRotation;
    lastRotation = currentRotation;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        predPosXs[i] = posXs[i] + velXs[i] * h;
        predPosYs[i] = posYs[i] + velYs[i] * h;
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
        fixed16 force_x, force_y;

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

        velXs[i] += force_x * h;
        velYs[i] += force_y * h;

        velXs[i] *= fixed16(static_cast<float>(1. - DAMPING * H));
        velYs[i] *= fixed16(static_cast<float>(1. - DAMPING * H));

        posXs[i] += velXs[i] * h;
        posYs[i] += velYs[i] * h;
    }

    // HACK

    fixed16 c = fixed16(cos(currentRotation.to_float()));
    fixed16 s = fixed16(sin(currentRotation.to_float()));

    fixed16 cd = fixed16(cos(dRotation.to_float()));
    fixed16 sd = fixed16(sin(dRotation.to_float()));

    fixed16 normalXs[] = { c, s, -c, -s };
    fixed16 normalYs[] = { s, -c, -s, c };
    fixed16 distances[] = {
        fixed16(0.5f * WIDTH),
        fixed16(0.5f * HEIGHT),
        fixed16(0.5f * WIDTH),
        fixed16(0.5f * HEIGHT)
    };

    // return;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        for (int w = 0; w < 4; w++)
        {
            fixed16 normalX = normalXs[w];
            fixed16 normalY = normalYs[w];

            fixed16 distance = distances[w];
            fixed16 pdot =
                normalX * (posXs[i] + normalX * distance) +
                normalY * (posYs[i] + normalY * distance);

            if (pdot < fixed16(0.f))
            {
                fixed16 wallVelX =
                    ((cd - fixed16(1.f)) * posXs[i] - sd * posYs[i]) / h;
                fixed16 wallVelY =
                    (sd * posXs[i] + (cd - fixed16(1.f)) * posYs[i]) / h;

                fixed16 relVelX = velXs[i] - wallVelX;
                fixed16 relVelY = velYs[i] - wallVelY;

                fixed16 vdot = normalX * relVelX + normalY * relVelY;

                fixed16 reflectedVelX =
                    relVelX - fixed16(static_cast<float>(1. + WALL_RESTITUTION)) * vdot * normalX + wallVelX;
                fixed16 reflectedVelY =
                    relVelY - fixed16(static_cast<float>(1. + WALL_RESTITUTION)) * vdot * normalY + wallVelY;

                fixed16 dtFrac = pdot / vdot;
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
    fixed16 posXs[NUM_PARTICLES] = { fixed16(0.f) };
    fixed16 posYs[NUM_PARTICLES] = { fixed16(0.f) };

    // HACK
    // although this only run once so it's not much of a concern

    float counter = 0;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        float x = fmod(counter, 1.f);
        float y = (counter - x) * 0.03f;

        posXs[i] = fixed16((x - 0.5f) * WIDTH);
        posYs[i] = fixed16((y - 0.5f) * HEIGHT);
        counter += 0.03;
    }

    fixed16 predPosXs[NUM_PARTICLES] = { fixed16(0.f) };
    fixed16 predPosYs[NUM_PARTICLES] = { fixed16(0.f) };

    fixed16 velXs[NUM_PARTICLES] = { fixed16(0.f) };
    fixed16 velYs[NUM_PARTICLES] = { fixed16(0.f) };

    fixed16 densities[NUM_PARTICLES] = { fixed16(0.f) };
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

        bool buf[DISPLAY_WIDTH][DISPLAY_HEIGHT] = { false };

        fixed16 c = fixed16(cos(-currentRotation.to_float()));
        fixed16 s = fixed16(sin(-currentRotation.to_float()));

        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            fixed16 x = c * posXs[i] - s * posYs[i];
            fixed16 y = s * posXs[i] + c * posYs[i];

            int xi = (
                (fixed16(0.5f) + x / fixed16(static_cast<float>(WIDTH)))
                * fixed16(static_cast<float>(DISPLAY_WIDTH))).to_int();
            int yi = (
                (fixed16(0.5f) + y / fixed16(static_cast<float>(HEIGHT)))
                * fixed16(static_cast<float>(DISPLAY_HEIGHT))).to_int();

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
