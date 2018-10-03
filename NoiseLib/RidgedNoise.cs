using System;

namespace Noise
{
    public class RidgedNoise
    {
        public static float[] Calc1D(int width, float scale)
        {
            float[] values = new float[width];
            for (int i = 0; i < width; i++)
                values[i] = Generate(i * scale) * 128 + 128;
            return values;
        }

        public static float[,] Calc2D(int width, int height, float scale)
        {
            float[,] values = new float[width, height];
            for (int i = 0; i < width; i++)
                for (int j = 0; j < height; j++)
                    values[i, j] = Generate(i * scale, j * scale) * 128 + 128;
            return values;
        }

        public static float[,,] Calc3D(int width, int height, int length, float scale)
        {
            float[,,] values = new float[width, height, length];
            for (int i = 0; i < width; i++)
                for (int j = 0; j < height; j++)
                    for (int k = 0; k < length; k++)
                        values[i, j, k] = Generate(i * scale, j * scale, k * scale) * 128 + 128;
            return values;
        }

        public static float CalcPixel1D(int x, float scale)
        {
            return Generate(x * scale) * 128 + 128;
        }

        public static float CalcPixel2D(int x, int y, float scale)
        {
            return Generate(x * scale, y * scale) * 128 + 128;
        }

        public static float CalcPixel3D(int x, int y, int z, float scale)
        {
            return Generate(x * scale, y * scale, z * scale) * 128 + 128;
        }

        static RidgedNoise()
        {
        }

        internal static float Generate(float x)
        {
            return 1.0f - Math.Abs(Perlin.Generate(x));
        }

        internal static float Generate(float x, float y)
        {
            return 1.0f - Math.Abs(Perlin.Generate(x, y));
        }

        internal static float Generate(float x, float y, float z)
        {
            return 1.0f - Math.Abs(Perlin.Generate(x, y, z));
        }

        internal static float Generate(float x, float y, float z, float w)
        {
            return 1.0f - Math.Abs(Perlin.Generate(x, y, z, w));
        }
    }
}
