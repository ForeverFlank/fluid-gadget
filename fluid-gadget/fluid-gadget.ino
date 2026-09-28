#include <SPI.h>
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SH110X.h>

#define i2c_Address 0x3c

#define SCREEN_WIDTH 128
#define SCREEN_HEIGHT 64
#define OLED_RESET -1
Adafruit_SH1106G display = Adafruit_SH1106G(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire, OLED_RESET);


#define NUMFLAKES 10
#define XPOS 0
#define YPOS 1
#define DELTAY 2

uint8_t display_buffer[1024] = {0};
int i = 0;

void setup()
{
    Serial.begin(9600);

    delay(250);
    display.begin(i2c_Address, true);

    // display.display();
    // delay(100);

    display.clearDisplay();
}

void loop()
{
    display_buffer[i]++;
    i++;
    if (i == 1024) i = 0;

    display.drawBitmap(0, 0, display_buffer, 128, 64, 1);
    display.display();

    delay(100);
}
