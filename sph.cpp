#include <iostream>
#include <unistd.h>
#include <sstream>
#include <cstdlib>
#include <cmath>

#define NUM_PARTICLES   64
#define WIDTH           16
#define HEIGHT          16

#define DT              0.01
#define KERNEL_RADIUS   5
#define PRESSURE_FORCE  5

using namespace std;

struct Vector2
{
    float x;
    float y;

    Vector2()
    {
        this->x = 0;
        this->y = 0;
    }

    Vector2(float x, float y)
    {
        this->x = x;
        this->y = y;
    }
};

float smoothingKernel(float distance)
{
    float normalized = distance / KERNEL_RADIUS;
    float value = max(1.f - normalized * normalized, 0.f);
    return value * value * value;
}

float calculateDensity(Vector2 positions[], float x, float y)
{
    float res = 0;

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        float dx = positions[i].x - x;
        float dy = positions[i].y - y;
        float distance = sqrt(dx * dx + dy * dy);
        res += smoothingKernel(distance);
    }

    return res;
}

float getDensityBounded(float density[][HEIGHT], int x, int y)
{
    x = min(max(x, 0), WIDTH - 1);
    y = min(max(y, 0), HEIGHT - 1);

    return density[x][y];
}

Vector2 calculatePressure(float density[][HEIGHT], int x, int y)
{
    float u = getDensityBounded(density, x, y - 1);
    float d = getDensityBounded(density, x, y + 1);
    float l = getDensityBounded(density, x - 1, y);
    float r = getDensityBounded(density, x + 1, y);

    return Vector2((r - l) * 0.5f, (d - u) * 0.5f);
}

int main()
{
    Vector2 positions[NUM_PARTICLES];
    Vector2 velocities[NUM_PARTICLES];
    float density[WIDTH][HEIGHT];

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        positions[i] = Vector2(i % WIDTH, i / WIDTH);
        velocities[i] = Vector2(0, 0);
    }

    for (int t = 0; t < 1000; t++)
    {
        for (int y = 0; y < HEIGHT; y++)
        {
            for (int x = 0; x < WIDTH; x++)
            {
                density[x][y] = calculateDensity(positions, x, y);
            }
        }

        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            float x = positions[i].x;
            float y = positions[i].y;
            Vector2 pressure = calculatePressure(density, (int)x, (int)y);

            velocities[i].x -= pressure.x * PRESSURE_FORCE * DT;
            velocities[i].y -= pressure.y * PRESSURE_FORCE * DT + 0.01f;

            positions[i].x += velocities[i].x * DT;
            positions[i].y += velocities[i].y * DT;

            positions[i].x = min(max(positions[i].x, 0.f), (float)WIDTH);
            positions[i].y = min(max(positions[i].y, 0.f), (float)HEIGHT);
        }

        if (t % 10 != 0) continue;

        char buf[WIDTH + 1][HEIGHT + 1] = {' '};

        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            int x = (int)positions[i].x;
            int y = HEIGHT - (int)positions[i].y - 1;

            buf[x][y] = ':';
        }

        stringstream ss;

        for (int y = 0; y < HEIGHT + 1; y++)
        {
            for (int x = 0; x < WIDTH + 1; x++)
            {
                ss << buf[x][y] << buf[x][y];
                // char ch[] = {' ', '.', ':', '#'};
                // ss << ch[min(4, (int)(density[x][y] * 5))];
                // ss << density[x][y] << " ";
                // ss << density[x][y] << " ";
            }
            ss << "\n";
        }

        usleep(10000);
        cout << "\033[2J\033[1;1H";
        cout << ss.str();
        cout << "----------\n";

        // for (int i = 0; i < NUM_PARTICLES; i++)
        // {
        //     cout << positions[i].x << " " << positions[i].y << "\n";
        // }
    }

    return 0;
}