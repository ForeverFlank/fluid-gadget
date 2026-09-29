#include <SPI.h>
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SH110X.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_MPU6050.h>
#include "../fixed16.cpp"


#define SCREEN_WIDTH        128
#define SCREEN_HEIGHT       64
#define OLED_RESET          -1
Adafruit_SH1106G display = Adafruit_SH1106G(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire, OLED_RESET);

Adafruit_MPU6050 mpu;


#define PIXEL_SIZE          4

#define WIDTH               50
#define HEIGHT              25

#define DISPLAY_WIDTH       SCREEN_WIDTH / PIXEL_SIZE
#define DISPLAY_HEIGHT      SCREEN_HEIGHT / PIXEL_SIZE

#define DT                  0.1
#define SUBSTEPS            4
#define H                   DT / SUBSTEPS

#define NUM_PARTICLES       64
#define PARTICLES_PER_CELL  4

#define KERNEL_RADIUS       1.25
#define PRESSURE_MULT       20
#define VISCOSITY_MULT      0.01
#define DAMPING             0.01
#define WALL_RESTITUTION    0.6

#define CONTAINER_HYP       sqrt(WIDTH * WIDTH + HEIGHT * HEIGHT)
#define SPATIAL_GRID_DIM    (int)ceil(CONTAINER_HYP / KERNEL_RADIUS)


const fixed16 dt = fixed16(static_cast<float>(DT));
const fixed16 h  = fixed16(static_cast<float>(H));

const fixed16 kernelRadius  = fixed16(static_cast<float>(KERNEL_RADIUS));
const fixed16 pressuerMult  = fixed16(static_cast<float>(PRESSURE_MULT));
const fixed16 viscosityMult = fixed16(static_cast<float>(VISCOSITY_MULT));
const fixed16 damping       = fixed16(static_cast<float>(DAMPING));
fixed16 reflectFactor = fixed16(static_cast<float>(1.f + WALL_RESTITUTION));


fixed16 gravityX = fixed16(0.f);
fixed16 gravityY = fixed16(-1.f);
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
    if (xi >= SPATIAL_GRID_DIM) xi = SPATIAL_GRID_DIM - 1;

    if (yi < 0) yi = 0;
    if (yi >= SPATIAL_GRID_DIM) yi = SPATIAL_GRID_DIM - 1;

    *out_xi = xi;
    *out_yi = yi;
}

void updateSpatialGrid(
    fixed16 *predPosXs,
    fixed16 *predPosYs,
    SpatialGridCell *spatialGrid)
{
    for (int yi = 0; yi < SPATIAL_GRID_DIM; yi++)
    {
        for (int xi = 0; xi < SPATIAL_GRID_DIM; xi++)
        {
            spatialGrid[yi * SPATIAL_GRID_DIM + xi].size = 0;
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

        SpatialGridCell &cell = spatialGrid[yi * SPATIAL_GRID_DIM + xi];

        if (cell.size == PARTICLES_PER_CELL - 1) continue;

        cell.array[cell.size++] = i;
        // spatialGrid[index].size += 1;
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
                out_indices[(*out_indicesSize)++] = cell.array[i];
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
        forceY += pressure * diry * pressuerMult;

        // viscosity
        fixed16 dvx = velXs[i] - velXs[index];
        fixed16 dvy = velYs[i] - velYs[index];

        fixed16 viscosityInfluence = smoothingKernel(dist);

        forceX += dvx * viscosityInfluence * viscosityMult;
        forceY += dvy * viscosityInfluence * viscosityMult;
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

        velXs[i] -= velXs[i] * damping * h;
        velYs[i] -= velYs[i] * damping * h;

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

            if (pdot >= fixed16(0.f)) continue;

            fixed16 wallVelX =
                ((cd - fixed16(1.f)) * posXs[i] - sd * posYs[i]) / h;
            fixed16 wallVelY =
                (sd * posXs[i] + (cd - fixed16(1.f)) * posYs[i]) / h;

            fixed16 relVelX = velXs[i] - wallVelX;
            fixed16 relVelY = velYs[i] - wallVelY;

            fixed16 vdot = normalX * relVelX + normalY * relVelY;

            if (vdot == fixed16(0.f)) continue;

            fixed16 reflectedVelX = relVelX - reflectFactor * vdot * normalX + wallVelX;
            fixed16 reflectedVelY = relVelY - reflectFactor * vdot * normalY + wallVelY;

            fixed16 dtFrac = pdot / vdot;

            posXs[i] += (reflectedVelX - velXs[i]) * dtFrac;
            posYs[i] += (reflectedVelY - velYs[i]) * dtFrac;

            velXs[i] = reflectedVelX;
            velYs[i] = reflectedVelY;
        }
    }
};


fixed16 posXs[NUM_PARTICLES] = { fixed16(0.f) };
fixed16 posYs[NUM_PARTICLES] = { fixed16(0.f) };

fixed16 predPosXs[NUM_PARTICLES] = { fixed16(0.f) };
fixed16 predPosYs[NUM_PARTICLES] = { fixed16(0.f) };

fixed16 velXs[NUM_PARTICLES] = { fixed16(0.f) };
fixed16 velYs[NUM_PARTICLES] = { fixed16(0.f) };

fixed16 densities[NUM_PARTICLES] = { fixed16(0.f) };
SpatialGridCell spatialGrid[SPATIAL_GRID_DIM * SPATIAL_GRID_DIM];


// #define BITMAP_SIZE SCREEN_WIDTH * SCREEN_HEIGHT / 8

// uint8_t bitmap[BITMAP_SIZE] = { 0 };
int i = 0;

void setup()
{
    Serial.begin(9600);

    delay(250);
    display.begin(0x3c, true);

    while (!mpu.begin())
    {
        Serial.println("Failed to find MPU6050 chip");
        delay(1000);
    }

    mpu.setAccelerometerRange(MPU6050_RANGE_8_G);
    mpu.setGyroRange(MPU6050_RANGE_500_DEG);
    mpu.setFilterBandwidth(MPU6050_BAND_21_HZ);

    float counter = 0;
    const float step = 0.021f;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        float x = fmod(counter, 0.51f);
        float y = (counter - x) * step * 5.f;

        posXs[i] = fixed16((x - 0.5f) * WIDTH);
        posYs[i] = fixed16((y - 0.5f) * HEIGHT);
        counter += step;
    }

    display.clearDisplay();
    display.display();
}

int oldXs[NUM_PARTICLES];
int oldYs[NUM_PARTICLES];

void mpuRead(int16_t *accel)
{
    uint8_t data[6] = { 0, 0, 0, 0, 0, 0 };
    Wire.beginTransmission(0x68);
    Wire.write(0x3b);
    Wire.endTransmission(false);
    Wire.requestFrom(0x68, 6, true);
    for (uint8_t i = 0; i < 6; ++i)
    {
        data[i] = Wire.read();
    }
    accel[0] = (int16_t)(data[0] << 8 | data[1]);
    accel[1] = (int16_t)(data[2] << 8 | data[3]);
    accel[2] = (int16_t)(data[4] << 8 | data[5]);

    // Wire.beginTransmission(0x68);
    // Wire.write(0x43);
    // Wire.endTransmission(false);
    // Wire.requestFrom(0x68, 6, true);
    // for (uint8_t i = 0; i < 6; ++i) {
    //     data[i] = Wire.read();
    // }
    // mpuArr[3] = (int16_t)(data[0] << 8 | data[1]);
    // mpuArr[4] = (int16_t)(data[2] << 8 | data[3]);
    // mpuArr[5] = (int16_t)(data[4] << 8 | data[5]);
    // return true;
}

bool firstFrame = true;
fixed16 sensorScale = fixed16(0.1f);

fixed16 prevAccelX;
fixed16 prevAccelY;
fixed16 prevAccelZ;

void loop()
{
    long startMillis = millis();

    sensors_event_t a, g, temp;
    mpu.getEvent(&a, &g, &temp);

    // HACK
    fixed16 accelX = sensorScale * fixed16(-a.acceleration.y);
    fixed16 accelY = sensorScale * fixed16(a.acceleration.x);
    fixed16 accelZ = sensorScale * fixed16(a.acceleration.z);

    gravityX = accelX;
    gravityY = accelY;

    if (!firstFrame)
    {
        fixed16 impulseX = (accelX - prevAccelX) / dt;
        fixed16 impulseY = (accelY - prevAccelY) / dt;

        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            velXs[i] += impulseX;
            velYs[i] += impulseY;
        }
    }

    prevAccelX = accelX;
    prevAccelY = accelY;
    prevAccelZ = accelZ;

    // Serial.print("Acceleration X: ");
    // Serial.print(a.acceleration.x);
    // Serial.print(", Y: ");
    // Serial.print(a.acceleration.y);
    // Serial.print(", Z: ");
    // Serial.print(a.acceleration.z);
    // Serial.println(" m/s^2");

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

    long updateSimMillis = millis();

    fixed16 c = fixed16(cos(-currentRotation.to_float()));
    fixed16 s = fixed16(sin(-currentRotation.to_float()));

    // for (int i = 0; i < BITMAP_SIZE; i++)
    // {
    //     bitmap[i] = 0;
    // }


    if (!firstFrame)
    {
        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            display.fillRect(oldXs[i] * PIXEL_SIZE, oldYs[i] * PIXEL_SIZE, PIXEL_SIZE, PIXEL_SIZE, 0);
        }
    }

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

        yi = DISPLAY_HEIGHT - yi - 1;

        // int bitmapIndex = (yi * DISPLAY_WIDTH + xi) / 8;

        // bitmap[bitmapIndex] |= 1 << (xi % 8);
        display.fillRect(xi * PIXEL_SIZE, yi * PIXEL_SIZE, PIXEL_SIZE, PIXEL_SIZE, 1);

        oldXs[i] = xi;
        oldYs[i] = yi;
    }

    display.display();


    long displayMillis = millis();


    // bitmap[i]++;
    // i++;
    // if (i == 1024) i = 0;

    // display.drawBitmap(0, 0, bitmap, DISPLAY_WIDTH, DISPLAY_HEIGHT, 1);

    Serial.print("update ms = ");
    Serial.print(updateSimMillis - startMillis);

    Serial.print("    display ms = ");
    Serial.println(displayMillis - updateSimMillis);

    firstFrame = false;
}
