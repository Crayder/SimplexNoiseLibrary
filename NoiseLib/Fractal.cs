using System;

namespace Noise
{
    public class Fractal
    {
        public static float Generate(float x, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;

            for (int i = 0; i < octaves; i++)
            {
                float n = Perlin.Generate(x * freq);
                sum += n * amp;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
        public static float Generate(float x, float y, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;

            for (int i = 0; i < octaves; i++)
            {
                float n = Perlin.Generate(x * freq, y * freq);
                sum += n * amp;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
        public static float Generate(float x, float y, float z, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f, float frequency = 1.0f)
        {
            float sum = 0.0f;
            float freq = frequency;
            float amp = 0.5f;

            for (int i = 0; i < octaves; i++)
            {
                float n = Perlin.Generate(x * freq, y * freq, z * freq);
                sum += n * amp;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
        public static float Generate(float x, float y, float z, float w, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;

            for (int i = 0; i < octaves; i++)
            {
                float n = Perlin.Generate(x * freq, y * freq, z * freq, w * freq);
                sum += n * amp;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
        public static float GenerateWorley(float x, float y, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;

            for (int i = 0; i < octaves; i++)
            {
                float n = Worley.Generate(x * freq, y * freq);
                sum += n * amp;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
        public static float GenerateWorley(float x, float y, float z, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;

            for (int i = 0; i < octaves; i++)
            {
                float n = Worley.Generate(x * freq, y * freq, z * freq);
                sum += n * amp;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
        public static float GenerateWorleySmooth(float x, float y, float falloff, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;

            for (int i = 0; i < octaves; i++)
            {
                float n = Worley.GenerateSmooth(x * freq, y * freq, falloff);
                sum += n * amp;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
        public static float GenerateWorleySmooth(float x, float y, float z, float falloff, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;

            for (int i = 0; i < octaves; i++)
            {
                float n = Worley.GenerateSmooth(x * freq, y * freq, z * freq, falloff);
                sum += n * amp;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
        public static void GenerateDerivatives(float x, out float rx, out float ry, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float freq = 1.0f;
            float amp = 0.5f;

            rx = ry = 0.0f;

            for (int i = 0; i < octaves; i++)
            {
                float nx, ny;

                Perlin.GenerateDerivatives(x * freq, out nx, out ny);

                rx += nx * amp;
                ry += ny * amp;

                freq *= lacunarity;
                amp *= gain;
            }
        }
        public static void GenerateDerivatives(float x, float y, out float rx, out float ry, out float rz, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float freq = 1.0f;
            float amp = 0.5f;

            rx = ry = rz = 0.0f;

            for (int i = 0; i < octaves; i++)
            {
                float nx, ny, nz;

                Perlin.GenerateDerivatives(x * freq, y * freq, out nx, out ny, out nz);

                rx += nx * amp;
                ry += ny * amp;
                rz += nz * amp;

                freq *= lacunarity;
                amp *= gain;
            }
        }
        public static void GenerateDerivatives(float x, float y, float z, out float rx, out float ry, out float rz, out float rw, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float freq = 1.0f;
            float amp = 0.5f;

            rx = ry = rz = rw = 0.0f;

            for (int i = 0; i < octaves; i++)
            {
                float nx, ny, nz, nw;
                nx = ny = nz = nw = 0.0f;

                Perlin.GenerateDerivatives(x * freq, y * freq, z * freq, out nx, out ny, out nz, out nw);

                rx += nx * amp;
                ry += ny * amp;
                rz += nz * amp;
                rw += nw * amp;

                freq *= lacunarity;
                amp *= gain;
            }
        }

        public static void GenerateCurl(float x, float y, out float rx, out float ry, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float dx, dy, dz;

            GenerateDerivatives(x, y, out dx, out dy, out dz, octaves, lacunarity, gain);
            rx = dz;
            ry = -dy;
        }

        private static float ridge(float h, float offset)
        {
            h = offset - Math.Abs(h);
            return h * h;
        }

        public static float GenerateRidged(float x, float ridgeOffset, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;
            float prev = 1.0f;

            for (int i = 0; i < octaves; i++)
            {
                float n = ridge(Perlin.Generate(x * freq), ridgeOffset);
                sum += n * amp * prev;
                prev = n;
                freq *= lacunarity;
                amp *= gain;
            }
            return sum;
        }

        public static float GenerateRidged(float x, float y, float ridgeOffset, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;
            float prev = 1.0f;

            for (int i = 0; i < octaves; i++)
            {
                float n = ridge(Perlin.Generate(x * freq, y * freq), ridgeOffset);
                sum += n * amp * prev;
                prev = n;
                freq *= lacunarity;
                amp *= gain;
            }
            return sum;
        }

        public static float GenerateRidged(float x, float y, float z, float ridgeOffset, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;
            float prev = 1.0f;

            for (int i = 0; i < octaves; i++)
            {
                float n = ridge(Perlin.Generate(x * freq, y * freq, z * freq), ridgeOffset);
                sum += n * amp * prev;
                prev = n;
                freq *= lacunarity;
                amp *= gain;
            }
            return sum;
        }

        public static float GenerateRidged(float x, float y, float z, float w, float ridgeOffset, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f)
        {
            float sum = 0.0f;
            float freq = 1.0f;
            float amp = 0.5f;
            float prev = 1.0f;

            for (int i = 0; i < octaves; i++)
            {
                float n = ridge(Perlin.Generate(x * freq, y * freq, z * freq, w * freq), ridgeOffset);
                sum += n * amp * prev;
                prev = n;
                freq *= lacunarity;
                amp *= gain;
            }

            return sum;
        }
    }
}
