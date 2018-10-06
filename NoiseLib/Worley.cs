using System;

namespace Noise
{
    public class Worley
    {
        public static float Generate(float x, float y)
        {
            float px = (float)Math.Floor(x), py = (float)Math.Floor(y);
            float fx = x - (float)Math.Floor(x), fy = y - (float)Math.Floor(y);

            float res = 8.0f;
            for (int j = -1; j <= 1; j++)
            {
                for (int i = -1; i <= 1; i++)
                {
                    float bx = i, by = j;
                    float n = Perlin.Generate(px + bx, py + by) * 0.5f + 0.5f;
                    float rx = bx - fx + n;
                    float ry = by - fy + n;
                    float d = (rx * rx + ry * ry);
                    res = Math.Min(res, d);
                }
            }
            return (float)Math.Sqrt(res);
        }

        public static float Generate(float x, float y, float z)
        {
            float px = (float)Math.Floor(x), py = (float)Math.Floor(y), pz = (float)Math.Floor(z);
            float fx = x - (float)Math.Floor(x), fy = y - (float)Math.Floor(y), fz = z - (float)Math.Floor(z);

            float res = 8.0f;
            for (int k = -1; k <= 1; k++)
            {
                for (int j = -1; j <= 1; j++)
                {
                    for (int i = -1; i <= 1; i++)
                    {
                        float bx = i, by = j, bz = k;
                        float n = Perlin.Generate(px + bx, py + by, pz + bz) * 0.5f + 0.5f;
                        float rx = bx - fx + n;
                        float ry = by - fy + n;
                        float rz = bz - fz + n;
                        float d = (rx * rx + ry * ry + rz * rz);
                        res = Math.Min(res, d);
                    }
                }
            }
            return (float)Math.Sqrt(res);
        }

        public static float GenerateSmooth(float x, float y, float falloff)
        {
            float px = (float)Math.Floor(x), py = (float)Math.Floor(y);
            float fx = x - (float)Math.Floor(x), fy = y - (float)Math.Floor(y);

            float res = 0.0f;
            for (int j = -1; j <= 1; j++)
            {
                for (int i = -1; i <= 1; i++)
                {
                    float bx = i, by = j;
                    float n = Perlin.Generate(px + bx, py + by) * 0.5f + 0.5f;
                    float rx = bx - fx + n;
                    float ry = by - fy + n;
                    float d = (rx * rx + ry * ry);
                    res += (float)Math.Exp(-falloff * d);
                }
            }
            return -(1.0f / falloff) * (float)Math.Log10(res);
        }

        public static float GenerateSmooth(float x, float y, float z, float falloff)
        {
            float px = (float)Math.Floor(x), py = (float)Math.Floor(y), pz = (float)Math.Floor(z);
            float fx = x - (float)Math.Floor(x), fy = y - (float)Math.Floor(y), fz = z - (float)Math.Floor(z);

            float res = 0.0f;
            for (int k = -1; k <= 1; k++)
            {
                for (int j = -1; j <= 1; j++)
                {
                    for (int i = -1; i <= 1; i++)
                    {
                        float bx = i, by = j, bz = k;
                        float n = Perlin.Generate(px + bx, py + by, pz + bz) * 0.5f + 0.5f;
                        float rx = bx - fx + n;
                        float ry = by - fy + n;
                        float rz = bz - fz + n;
                        float d = (rx * rx + ry * ry + rz * rz);
                        res += (float)Math.Exp(-falloff * d);
                    }
                }
            }
            return -(1.0f / falloff) * (float)Math.Log10(res);
        }
    }
}
