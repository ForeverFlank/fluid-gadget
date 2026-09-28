#include <SPI.h>
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SH110X.h>

#define i2c_Address 0x3c

#define SCREEN_WIDTH        128
#define SCREEN_HEIGHT       64
#define OLED_RESET          -1
Adafruit_SH1106G display = Adafruit_SH1106G(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire, OLED_RESET);


#define WIDTH               128
#define HEIGHT              64

#define DT                  0.1
#define SUBSTEPS            1
#define H                   DT / SUBSTEPS

#define NUM_PARTICLES       256
#define PARTICLES_PER_CELL  8

#define KERNEL_RADIUS       16
#define WALL_RESTITUTION    0.85
#define PRESSURE_MULT       20
#define VISCOSITY_MULT      0
#define DAMPING             0.0

#define CONTAINER_HYP       sqrt(WIDTH * WIDTH + HEIGHT * HEIGHT)
#define SPATIAL_GRID_DIM    (int)ceil(CONTAINER_HYP / KERNEL_RADIUS)


float gravityX = 0;
float gravityY = 10;
float currentRotation = 0;
float lastRotation = currentRotation;


struct SpatialGridCell
{
    int array[PARTICLES_PER_CELL];
    uint16_t size;
};


float smoothingKernel(float distance)
{
    float normalized = distance / KERNEL_RADIUS;
    float value = max(1.f - normalized * normalized, 0.f);
    return value * value * value;
}

float d_smoothingKernel(float distance)
{
    float normalized = distance / KERNEL_RADIUS;
    float value = max(1 - normalized * normalized, 0.f);
    return -6.f * value * value * normalized / KERNEL_RADIUS;
}


void spatialGridCoord(
    float *predPosXs,
    float *predPosYs,
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
    if (xi >= SPATIAL_GRID_DIM) xi = SPATIAL_GRID_DIM - 1;

    if (yi < 0) yi = 0;
    if (yi >= SPATIAL_GRID_DIM) yi = SPATIAL_GRID_DIM - 1;

    *out_xi = xi;
    *out_yi = yi;
}

void updateSpatialGrid(
    float *predPosXs,
    float *predPosYs,
    SpatialGridCell *spatialGrid)
{
    for (int i = 0; i < SPATIAL_GRID_DIM; i++)
    {
        for (int j = 0; j < SPATIAL_GRID_DIM; j++)
        {
            spatialGrid[j * SPATIAL_GRID_DIM + i].size = 0;
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

        if (spatialGrid[yi * SPATIAL_GRID_DIM + xi].size == PARTICLES_PER_CELL - 1) continue;

        spatialGrid[yi * SPATIAL_GRID_DIM + xi].array[
            spatialGrid[yi * SPATIAL_GRID_DIM + xi].size
        ] = i;

        spatialGrid[yi * SPATIAL_GRID_DIM + xi].size += 1;
    }
};

void nearbyIndices(
    float *predPosXs,
    float *predPosYs,
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
    int xiMax = min(xi + 1, SPATIAL_GRID_DIM - 1);
    int yiMin = max(yi - 1, 0);
    int yiMax = min(yi + 1, SPATIAL_GRID_DIM - 1);

    for (int xj = xiMin; xj <= xiMax; xj++)
    {
        for (int yj = yiMin; yj <= yiMax; yj++)
        {
            SpatialGridCell cell = spatialGrid[yj * SPATIAL_GRID_DIM + xj];

            for (int i = 0; i < cell.size; i++)
            {
                out_indices[*out_indicesSize] = cell.array[i];
                (*out_indicesSize)++;
            }
        }
    }
}


float calculateDensity(
    float *predPosXs,
    float *predPosYs,
    SpatialGridCell *spatialGrid,
    int index)
{
    float res = 0;

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

        float dx = predPosXs[i] - predPosXs[index];
        float dy = predPosYs[i] - predPosYs[index];
        float dist = sqrt(dx * dx + dy * dy);

        if (dist > KERNEL_RADIUS) continue;

        res += smoothingKernel(dist);
    }

    return res;
}


void calculateForce(
    float *predPosXs,
    float *predPosYs,
    float *velXs,
    float *velYs,
    float *densities,
    SpatialGridCell *spatialGrid,
    int index,
    float *out_forceX,
    float *out_forceY)
{
    float forceX = gravityX;
    float forceY = gravityY;

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

        float dx = predPosXs[i] - predPosXs[index];
        float dy = predPosYs[i] - predPosYs[index];
        float dist = sqrt(dx * dx + dy * dy);

        if (dist > KERNEL_RADIUS) continue;

        float dirx = dist > 0 ? dx / dist : 0;
        float diry = dist > 0 ? dy / dist : 0;

        float density = densities[i];
        float pressureSlope = d_smoothingKernel(dist);
        float pressure = density * pressureSlope;

        forceX += pressure * dirx * PRESSURE_MULT;
        forceY += pressure * diry * PRESSURE_MULT;

        // viscosity
        float dvx = velXs[i] - velXs[index];
        float dvy = velYs[i] - velYs[index];

        float viscosityInfluence = smoothingKernel(dist);

        forceX += dvx * viscosityInfluence * VISCOSITY_MULT;
        forceY += dvy * viscosityInfluence * VISCOSITY_MULT;
    }

    *out_forceX = forceX;
    *out_forceY = forceY;
}

void updateSim(
    float *posXs,
    float *posYs,
    float *predPosXs,
    float *predPosYs,
    float *velXs,
    float *velYs,
    float *densities,
    SpatialGridCell *spatialGrid)
{
    float dRotation = currentRotation - lastRotation;
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
        float force_x, force_y;

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


    float c = cos(currentRotation);
    float s = sin(currentRotation);

    float cd = cos(dRotation);
    float sd = sin(dRotation);

    float normalXs[] = {c, s, -c, -s};
    float normalYs[] = {s, -c, -s, c};
    float distances[] = {0.5f * WIDTH, 0.5f * HEIGHT, 0.5f * WIDTH, 0.5f * HEIGHT};

    // return;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        for (int w = 0; w < 4; w++)
        {
            float normalX = normalXs[w];
            float normalY = normalYs[w];

            float distance = distances[w];
            float pdot =
                normalX * (posXs[i] + normalX * distance) +
                normalY * (posYs[i] + normalY * distance);

            if (pdot < 0)
            {
                float wallVelX = ((cd - 1) * posXs[i] - sd * posYs[i]) / H;
                float wallVelY = (sd * posXs[i] + (cd - 1) * posYs[i]) / H;

                float relVelX = velXs[i] - wallVelX;
                float relVelY = velYs[i] - wallVelY;

                float vdot = normalX * relVelX + normalY * relVelY;

                float reflectedVelX =
                    relVelX - (1 + WALL_RESTITUTION) * vdot * normalX + wallVelX;
                float reflectedVelY =
                    relVelY - (1 + WALL_RESTITUTION) * vdot * normalY + wallVelY;

                float dtFrac = pdot / vdot;
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


float posXs[NUM_PARTICLES] = {0};
float posYs[NUM_PARTICLES] = {0};

float counter = 0;

float predPosXs[NUM_PARTICLES] = {0};
float predPosYs[NUM_PARTICLES] = {0};

float velXs[NUM_PARTICLES] = {0};
float velYs[NUM_PARTICLES] = {0};

float densities[NUM_PARTICLES] = {0};
SpatialGridCell spatialGrid[SPATIAL_GRID_DIM * SPATIAL_GRID_DIM];

#define BITMAP_SIZE         WIDTH * HEIGHT / 8

uint8_t bitmap[BITMAP_SIZE] = {0};
int i = 0;

void setup()
{
    Serial.begin(9600);

    delay(250);
    display.begin(i2c_Address, true);

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        float x = fmod(counter, 1.f);
        float y = (counter - x) * 0.03f + 0.2f;

        posXs[i] = (x - 0.5f) * WIDTH;
        posYs[i] = (y - 0.5f) * HEIGHT;
        counter += 0.03;
    }
}

void loop()
{
    for (int i = 0; i < SUBSTEPS; i++)
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
    }

    float c = cos(-currentRotation);
    float s = sin(-currentRotation);

    for (int i = 0; i < BITMAP_SIZE; i++)
    {
        bitmap[i] = 0;
    }

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        float x = c * posXs[i] - s * posYs[i];
        float y = s * posXs[i] + c * posYs[i];

        int xi = (int)((0.5f + x / WIDTH) * WIDTH);
        int yi = (int)((0.5f + y / HEIGHT) * HEIGHT);

        xi = min(xi, WIDTH - 1);
        yi = min(yi, HEIGHT - 1);

        int bitmapIndex = (yi * WIDTH + xi) / 8;

        bitmap[bitmapIndex] |= 1 << (xi % 8);
    }


    display.clearDisplay();

    // bitmap[i]++;
    // i++;
    // if (i == 1024) i = 0;

    display.drawBitmap(0, 0, bitmap, 128, 64, 1);
    display.display();

    // Serial.println(millis());
}
