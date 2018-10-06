using System;

namespace Noise
{
    public class Ridged
    {
        public static float[] Calc(int width, float scale = 1.0f)
        {
            float[] values = new float[width];
            for (int i = 0; i < width; i++)
                values[i] = Generate(i * scale);
            return values;
        }

        public static float[,] Calc(int width, int height, float scale = 1.0f)
        {
            float[,] values = new float[width, height];
            for (int i = 0; i < width; i++)
                for (int j = 0; j < height; j++)
                    values[i, j] = Generate(i * scale, j * scale);
            return values;
        }

        public static float[,,] Calc(int width, int height, int length, float scale = 1.0f)
        {
            float[,,] values = new float[width, height, length];
            for (int i = 0; i < width; i++)
                for (int j = 0; j < height; j++)
                    for (int k = 0; k < length; k++)
                        values[i, j, k] = Generate(i * scale, j * scale, k * scale);
            return values;
        }

        public static float CalcPixel(int x, float scale = 1.0f)
        {
            return Generate(x * scale);
        }

        public static float CalcPixel(int x, int y, float scale = 1.0f)
        {
            return Generate(x * scale, y * scale);
        }

        public static float CalcPixel(int x, int y, int z, float scale = 1.0f)
        {
            return Generate(x * scale, y * scale, z * scale);
        }

        static Ridged()
        {
        }

        public static float Generate(float x)
        {
            return 1.0f - Math.Abs(Perlin.Generate(x));
        }

        public static float Generate(float x, float y)
        {
            return 1.0f - Math.Abs(Perlin.Generate(x, y));
        }

        public static float Generate(float x, float y, float z)
        {
            return 1.0f - Math.Abs(Perlin.Generate(x, y, z));
        }

        public static float Generate(float x, float y, float z, float w)
        {
            return 1.0f - Math.Abs(Perlin.Generate(x, y, z, w));
        }
    }
}
