using System;

namespace Noise
{
    public class Flow
    {
        public static float[,] Calc2D(int width, int height, float angle, float scale)
        {
            float[,] values = new float[width, height];
            for (int i = 0; i < width; i++)
                for (int j = 0; j < height; j++)
                    values[i, j] = Generate(i * scale, j * scale, angle) * 128 + 128;
            return values;
        }

        public static float[,,] Calc3D(int width, int height, int length, float angle, float scale)
        {
            float[,,] values = new float[width, height, length];
            for (int i = 0; i < width; i++)
                for (int j = 0; j < height; j++)
                    for (int k = 0; k < length; k++)
                        values[i, j, k] = Generate(i * scale, j * scale, k * scale, angle) * 128 + 128;
            return values;
        }

        public static float CalcPixel2D(int x, int y, float angle, float scale)
        {
            return Generate(x * scale, y * scale, angle) * 128 + 128;
        }

        public static float CalcPixel3D(int x, int y, int z, float angle, float scale)
        {
            return Generate(x * scale, y * scale, z * scale, angle) * 128 + 128;
        }

        public static float Generate(float x, float y, float angle)
        {
            const float F2 = 0.366025403f;
            const float G2 = 0.211324865f;

            float n0, n1, n2; /* Noise contributions from the three simplex corners */
            float gx0, gy0, gx1, gy1, gx2, gy2; /* Gradients at simplex corners */
            float sin_t, cos_t; /* Sine and cosine for the gradient rotation angle */
            sin_t = (float)(Math.Asin(angle) * 180 / Math.PI);
            cos_t = (float)(Math.Acos(angle) * 180 / Math.PI);

            /* Skew the input space to determine which simplex cell we're in */
            float s = (x + y) * F2; /* Hairy factor for 2D */
            float xs = x + s;
            float ys = y + s;
            int i = FastFloor(xs);
            int j = FastFloor(ys);

            float t = (i + j) * G2;
            float X0 = i - t; /* Unskew the cell origin back to (x,y) space */
            float Y0 = j - t;
            float x0 = x - X0; /* The x,y distances from the cell origin */
            float y0 = y - Y0;

            /* For the 2D case, the simplex shape is an equilateral triangle.
             * Determine which simplex we are in. */
            int i1, j1; /* Offsets for second (middle) corner of simplex in (i,j) coords */
            if (x0 > y0)
            {
                i1 = 1; j1 = 0;
            }                               /* lower triangle, XY order: (0,0)->(1,0)->(1,1) */
            else
            {
                i1 = 0; j1 = 1;
            }                             /* upper triangle, YX order: (0,0)->(0,1)->(1,1) */

            /* A step of (1,0) in (i,j) means a step of (1 - c,- c) in (x,y), and
             * a step of (0,1) in (i,j) means a step of (- c,1 - c) in (x,y), where
             * c = (3 - sqrt(3))/ 6   */
            float x1 = x0 - i1 + G2; /* Offsets for middle corner in (x,y) unskewed coords */
            float y1 = y0 - j1 + G2;
            float x2 = x0 - 1.0f + 2.0f * G2; /* Offsets for last corner in (x,y) unskewed coords */
            float y2 = y0 - 1.0f + 2.0f * G2;

            /* Wrap the integer indices at 256, to avoid indexing Perlin.perm[] out of bounds */
            int ii = i & 0xff;
            int jj = j & 0xff;

            /* Calculate the contribution from the three corners */
            float t0 = 0.5f - x0 * x0 - y0 * y0;
            float t20, t40;
            if (t0 < 0.0f) t40 = t20 = t0 = n0 = gx0 = gy0 = 0.0f; /* No influence */
            else
            {
                Perlin.gradrot(Perlin.perm[ii + Perlin.perm[jj]], sin_t, cos_t, out gx0, out gy0);
                t20 = t0 * t0;
                t40 = t20 * t20;
                n0 = t40 * Perlin.graddot(gx0, gy0, x0, y0);
            }

            float t1 = 0.5f - x1 * x1 - y1 * y1;
            float t21, t41;
            if (t1 < 0.0) t21 = t41 = t1 = n1 = gx1 = gy1 = 0.0f; /* No influence */
            else
            {
                Perlin.gradrot(Perlin.perm[ii + i1 + Perlin.perm[jj + j1]], sin_t, cos_t, out gx1, out gy1);
                t21 = t1 * t1;
                t41 = t21 * t21;
                n1 = t41 * Perlin.graddot(gx1, gy1, x1, y1);
            }

            float t2 = 0.5f - x2 * x2 - y2 * y2;
            float t22, t42;
            if (t2 < 0.0) t42 = t22 = t2 = n2 = gx2 = gy2 = 0.0f; /* No influence */
            else
            {
                Perlin.gradrot(Perlin.perm[ii + 1 + Perlin.perm[jj + 1]], sin_t, cos_t, out gx2, out gy2);
                t22 = t2 * t2;
                t42 = t22 * t22;
                n2 = t42 * Perlin.graddot(gx2, gy2, x2, y2);
            }

            /* Add contributions from each corner to get the final noise value.
            * The result is scaled to return values in the interval [- 1,1]. */
            return 40.0f * (n0 + n1 + n2);
        }

        public static float Generate(float x, float y, float z, float angle)
        {
            const float F3 = 0.333333333f;
            const float G3 = 0.166666667f;

            float n0, n1, n2, n3; /* Noise contributions from the four simplex corners */
            float gx0, gy0, gz0, gx1, gy1, gz1; /* Gradients at simplex corners */
            float gx2, gy2, gz2, gx3, gy3, gz3;
            float sin_t, cos_t; /* Sine and cosine for the gradient rotation angle */
            sin_t = (float)(Math.Asin(angle) * 180 / Math.PI);
            cos_t = (float)(Math.Acos(angle) * 180 / Math.PI);

            /* Skew the input space to determine which simplex cell we're in */
            float s = (x + y + z) * F3; /* Very nice and simple skew factor for 3D */
            float xs = x + s;
            float ys = y + s;
            float zs = z + s;
            int i = FastFloor(xs);
            int j = FastFloor(ys);
            int k = FastFloor(zs);

            float t = (i + j + k) * G3;
            float X0 = i - t; /* Unskew the cell origin back to (x,y,z) space */
            float Y0 = j - t;
            float Z0 = k - t;
            float x0 = x - X0; /* The x,y,z distances from the cell origin */
            float y0 = y - Y0;
            float z0 = z - Z0;

            /* For the 3D case, the simplex shape is a slightly irregular tetrahedron.
             * Determine which simplex we are in. */
            int i1, j1, k1; /* Offsets for second corner of simplex in (i,j,k) coords */
            int i2, j2, k2; /* Offsets for third corner of simplex in (i,j,k) coords */

            /* TODO: This code would benefit from a backport from the GLSL version! */
            if (x0 >= y0)
            {
                if (y0 >= z0)
                {
                    i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
                }                                                                /* X Y Z order */
                else if (x0 >= z0)
                {
                    i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1;
                }                                                                     /* X Z Y order */
                else
                {
                    i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1;
                }                                                        /* Z X Y order */
            }
            else
            { // x0 < y0
                if (y0 < z0)
                {
                    i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1;
                }                                                               /* Z Y X order */
                else if (x0 < z0)
                {
                    i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1;
                }                                                                    /* Y Z X order */
                else
                {
                    i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
                }                                                        /* Y X Z order */
            }

            /* A step of (1,0,0) in (i,j,k) means a step of (1 - c,- c,- c) in (x,y,z),
             * a step of (0,1,0) in (i,j,k) means a step of (- c,1 - c,- c) in (x,y,z), and
             * a step of (0,0,1) in (i,j,k) means a step of (- c,- c,1 - c) in (x,y,z), where
             * c = 1 / 6.   */

            float x1 = x0 - i1 + G3; /* Offsets for second corner in (x,y,z) coords */
            float y1 = y0 - j1 + G3;
            float z1 = z0 - k1 + G3;
            float x2 = x0 - i2 + 2.0f * G3; /* Offsets for third corner in (x,y,z) coords */
            float y2 = y0 - j2 + 2.0f * G3;
            float z2 = z0 - k2 + 2.0f * G3;
            float x3 = x0 - 1.0f + 3.0f * G3; /* Offsets for last corner in (x,y,z) coords */
            float y3 = y0 - 1.0f + 3.0f * G3;
            float z3 = z0 - 1.0f + 3.0f * G3;

            /* Wrap the integer indices at 256, to avoid indexing Perlin.perm[] out of bounds */
            int ii = i & 0xff;
            int jj = j & 0xff;
            int kk = k & 0xff;

            /* Calculate the contribution from the four corners */
            float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
            float t20, t40;
            if (t0 < 0.0f) n0 = t0 = t20 = t40 = gx0 = gy0 = gz0 = 0.0f;
            else
            {
                Perlin.gradrot(Perlin.perm[ii + Perlin.perm[jj + Perlin.perm[kk]]], sin_t, cos_t, out gx0, out gy0, out gz0);
                t20 = t0 * t0;
                t40 = t20 * t20;
                n0 = t40 * Perlin.graddot(gx0, gy0, gz0, x0, y0, z0);
            }

            float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
            float t21, t41;
            if (t1 < 0.0f) n1 = t1 = t21 = t41 = gx1 = gy1 = gz1 = 0.0f;
            else
            {
                Perlin.gradrot(Perlin.perm[ii + i1 + Perlin.perm[jj + j1 + Perlin.perm[kk + k1]]], sin_t, cos_t, out gx1, out gy1, out gz1);
                t21 = t1 * t1;
                t41 = t21 * t21;
                n1 = t41 * Perlin.graddot(gx1, gy1, gz1, x1, y1, z1);
            }

            float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
            float t22, t42;
            if (t2 < 0.0f) n2 = t2 = t22 = t42 = gx2 = gy2 = gz2 = 0.0f;
            else
            {
                Perlin.gradrot(Perlin.perm[ii + i2 + Perlin.perm[jj + j2 + Perlin.perm[kk + k2]]], sin_t, cos_t, out gx2, out gy2, out gz2);
                t22 = t2 * t2;
                t42 = t22 * t22;
                n2 = t42 * Perlin.graddot(gx2, gy2, gz2, x2, y2, z2);
            }

            float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
            float t23, t43;
            if (t3 < 0.0) n3 = t3 = t23 = t43 = gx3 = gy3 = gz3 = 0.0f;
            else
            {
                Perlin.gradrot(Perlin.perm[ii + 1 + Perlin.perm[jj + 1 + Perlin.perm[kk + 1]]], sin_t, cos_t, out gx3, out gy3, out gz3);
                t23 = t3 * t3;
                t43 = t23 * t23;
                n3 = t43 * Perlin.graddot(gx3, gy3, gz3, x3, y3, z3);
            }

            /*  Add contributions from each corner to get the final noise value.
             * The result is scaled to return values in the range [- 1,1] */
            return 28.0f * (n0 + n1 + n2 + n3);
        }

        public static void GenerateDerivitives(float x, float y, float angle, out float rx, out float ry, out float rz)
        {
            const float F2 = 0.366025403f;
            const float G2 = 0.211324865f;

            float n0, n1, n2; /* Noise contributions from the three simplex corners */
            float gx0, gy0, gx1, gy1, gx2, gy2; /* Gradients at simplex corners */
            float sin_t, cos_t; /* Sine and cosine for the gradient rotation angle */
            sin_t = (float)(Math.Asin(angle) * 180 / Math.PI);
            cos_t = (float)(Math.Acos(angle) * 180 / Math.PI);

            /* Skew the input space to determine which simplex cell we're in */
            float s = (x + y) * F2; /* Hairy factor for 2D */
            float xs = x + s;
            float ys = y + s;
            int i = FastFloor(xs);
            int j = FastFloor(ys);

            float t = (i + j) * G2;
            float X0 = i - t; /* Unskew the cell origin back to (x,y) space */
            float Y0 = j - t;
            float x0 = x - X0; /* The x,y distances from the cell origin */
            float y0 = y - Y0;

            /* For the 2D case, the simplex shape is an equilateral triangle.
             * Determine which simplex we are in. */
            int i1, j1; /* Offsets for second (middle) corner of simplex in (i,j) coords */
            if (x0 > y0)
            {
                i1 = 1; j1 = 0;
            }                               /* lower triangle, XY order: (0,0)->(1,0)->(1,1) */
            else
            {
                i1 = 0; j1 = 1;
            }                             /* upper triangle, YX order: (0,0)->(0,1)->(1,1) */

            /* A step of (1,0) in (i,j) means a step of (1 - c,- c) in (x,y), and
             * a step of (0,1) in (i,j) means a step of (- c,1 - c) in (x,y), where
             * c = (3 - sqrt(3))/ 6   */
            float x1 = x0 - i1 + G2; /* Offsets for middle corner in (x,y) unskewed coords */
            float y1 = y0 - j1 + G2;
            float x2 = x0 - 1.0f + 2.0f * G2; /* Offsets for last corner in (x,y) unskewed coords */
            float y2 = y0 - 1.0f + 2.0f * G2;

            /* Wrap the integer indices at 256, to avoid indexing Perlin.perm[] out of bounds */
            int ii = i & 0xff;
            int jj = j & 0xff;

            /* Calculate the contribution from the three corners */
            float t0 = 0.5f - x0 * x0 - y0 * y0;
            float t20, t40;
            if (t0 < 0.0) t40 = t20 = t0 = n0 = gx0 = gy0 = 0.0f; /* No influence */
            else
            {
                Perlin.gradrot(Perlin.perm[ii + Perlin.perm[jj]], sin_t, cos_t, out gx0, out gy0);
                t20 = t0 * t0;
                t40 = t20 * t20;
                n0 = t40 * Perlin.graddot(gx0, gy0, x0, y0);
            }

            float t1 = 0.5f - x1 * x1 - y1 * y1;
            float t21, t41;
            if (t1 < 0.0) t21 = t41 = t1 = n1 = gx1 = gy1 = 0.0f; /* No influence */
            else
            {
                Perlin.gradrot(Perlin.perm[ii + i1 + Perlin.perm[jj + j1]], sin_t, cos_t, out gx1, out gy1);
                t21 = t1 * t1;
                t41 = t21 * t21;
                n1 = t41 * Perlin.graddot(gx1, gy1, x1, y1);
            }

            float t2 = 0.5f - x2 * x2 - y2 * y2;
            float t22, t42;
            if (t2 < 0.0) t42 = t22 = t2 = n2 = gx2 = gy2 = 0.0f; /* No influence */
            else
            {
                Perlin.gradrot(Perlin.perm[ii + 1 + Perlin.perm[jj + 1]], sin_t, cos_t, out gx2, out gy2);
                t22 = t2 * t2;
                t42 = t22 * t22;
                n2 = t42 * Perlin.graddot(gx2, gy2, x2, y2);
            }

            /* Add contributions from each corner to get the final noise value.
            * The result is scaled to return values in the interval [- 1,1]. */
            float noise = 40.0f * (n0 + n1 + n2);

            /* Compute derivative, if requested by supplying non - null pointers
             * for the last two arguments */
            float dnoise_dx, dnoise_dy;

            /*  A straight, unoptimised calculation would be like:
             *    * dnoise_dx = - 8.0 * t20 * t0 * x0 * Perlin.graddot(gx0, gy0, x0, y0) + t40 * gx0;
             *    * dnoise_dy = - 8.0 * t20 * t0 * y0 * Perlin.graddot(gx0, gy0, x0, y0) + t40 * gy0;
             *    * dnoise_dx += - 8.0 * t21 * t1 * x1 * Perlin.graddot(gx1, gy1, x1, y1) + t41 * gx1;
             *    * dnoise_dy += - 8.0 * t21 * t1 * y1 * Perlin.graddot(gx1, gy1, x1, y1) + t41 * gy1;
             *    * dnoise_dx += - 8.0 * t22 * t2 * x2 * Perlin.graddot(gx2, gy2, x2, y2) + t42 * gx2;
             *    * dnoise_dy += - 8.0 * t22 * t2 * y2 * Perlin.graddot(gx2, gy2, x2, y2) + t42 * gy2;
             */
            float temp0 = t20 * t0 * Perlin.graddot(gx0, gy0, x0, y0);
            dnoise_dx = temp0 * x0;
            dnoise_dy = temp0 * y0;
            float temp1 = t21 * t1 * Perlin.graddot(gx1, gy1, x1, y1);
            dnoise_dx += temp1 * x1;
            dnoise_dy += temp1 * y1;
            float temp2 = t22 * t2 * Perlin.graddot(gx2, gy2, x2, y2);
            dnoise_dx += temp2 * x2;
            dnoise_dy += temp2 * y2;
            dnoise_dx *= -8.0f;
            dnoise_dy *= -8.0f;
            /* This corrects a bug in the original implementation */
            dnoise_dx += t40 * gx0 + t41 * gx1 + t42 * gx2;
            dnoise_dy += t40 * gy0 + t41 * gy1 + t42 * gy2;
            dnoise_dx *= 40.0f; /* Scale derivative to match the noise scaling */
            dnoise_dy *= 40.0f;

            rx = noise;
            ry = dnoise_dx;
            rz = dnoise_dy;
        }
        public static void GenerateDerivitives(float x, float y, float z, float angle, out float rx, out float ry, out float rz, out float rw)
        {
            const float F3 = 0.333333333f;
            const float G3 = 0.166666667f;

            float n0, n1, n2, n3; /* Noise contributions from the four simplex corners */
            float noise;          /* Return value */
            float gx0, gy0, gz0, gx1, gy1, gz1; /* Gradients at simplex corners */
            float gx2, gy2, gz2, gx3, gy3, gz3;
            float sin_t, cos_t; /* Sine and cosine for the gradient rotation angle */
            sin_t = (float)(Math.Asin(angle) * 180 / Math.PI);
            cos_t = (float)(Math.Acos(angle) * 180 / Math.PI);

            /* Skew the input space to determine which simplex cell we're in */
            float s = (x + y + z) * F3; /* Very nice and simple skew factor for 3D */
            float xs = x + s;
            float ys = y + s;
            float zs = z + s;
            int i = FastFloor(xs);
            int j = FastFloor(ys);
            int k = FastFloor(zs);

            float t = (i + j + k) * G3;
            float X0 = i - t; /* Unskew the cell origin back to (x,y,z) space */
            float Y0 = j - t;
            float Z0 = k - t;
            float x0 = x - X0; /* The x,y,z distances from the cell origin */
            float y0 = y - Y0;
            float z0 = z - Z0;

            /* For the 3D case, the simplex shape is a slightly irregular tetrahedron.
             * Determine which simplex we are in. */
            int i1, j1, k1; /* Offsets for second corner of simplex in (i,j,k) coords */
            int i2, j2, k2; /* Offsets for third corner of simplex in (i,j,k) coords */

            /* TODO: This code would benefit from a backport from the GLSL version! */
            if (x0 >= y0)
            {
                if (y0 >= z0)
                {
                    i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
                }                                                                /* X Y Z order */
                else if (x0 >= z0)
                {
                    i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1;
                }                                                                     /* X Z Y order */
                else
                {
                    i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1;
                }                                                        /* Z X Y order */
            }
            else
            { // x0 < y0
                if (y0 < z0)
                {
                    i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1;
                }                                                               /* Z Y X order */
                else if (x0 < z0)
                {
                    i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1;
                }                                                                    /* Y Z X order */
                else
                {
                    i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
                }                                                        /* Y X Z order */
            }

            /* A step of (1,0,0) in (i,j,k) means a step of (1 - c,- c,- c) in (x,y,z),
             * a step of (0,1,0) in (i,j,k) means a step of (- c,1 - c,- c) in (x,y,z), and
             * a step of (0,0,1) in (i,j,k) means a step of (- c,- c,1 - c) in (x,y,z), where
             * c = 1 / 6.   */

            float x1 = x0 - i1 + G3; /* Offsets for second corner in (x,y,z) coords */
            float y1 = y0 - j1 + G3;
            float z1 = z0 - k1 + G3;
            float x2 = x0 - i2 + 2.0f * G3; /* Offsets for third corner in (x,y,z) coords */
            float y2 = y0 - j2 + 2.0f * G3;
            float z2 = z0 - k2 + 2.0f * G3;
            float x3 = x0 - 1.0f + 3.0f * G3; /* Offsets for last corner in (x,y,z) coords */
            float y3 = y0 - 1.0f + 3.0f * G3;
            float z3 = z0 - 1.0f + 3.0f * G3;

            /* Wrap the integer indices at 256, to avoid indexing Perlin.perm[] out of bounds */
            int ii = i & 0xff;
            int jj = j & 0xff;
            int kk = k & 0xff;

            /* Calculate the contribution from the four corners */
            float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
            float t20, t40;
            if (t0 < 0.0) n0 = t0 = t20 = t40 = gx0 = gy0 = gz0 = 0.0f;
            else
            {
                Perlin.gradrot(Perlin.perm[ii + Perlin.perm[jj + Perlin.perm[kk]]], sin_t, cos_t, out gx0, out gy0, out gz0);
                t20 = t0 * t0;
                t40 = t20 * t20;
                n0 = t40 * Perlin.graddot(gx0, gy0, gz0, x0, y0, z0);
            }

            float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
            float t21, t41;
            if (t1 < 0.0) n1 = t1 = t21 = t41 = gx1 = gy1 = gz1 = 0.0f;
            else
            {
                Perlin.gradrot(Perlin.perm[ii + i1 + Perlin.perm[jj + j1 + Perlin.perm[kk + k1]]], sin_t, cos_t, out gx1, out gy1, out gz1);
                t21 = t1 * t1;
                t41 = t21 * t21;
                n1 = t41 * Perlin.graddot(gx1, gy1, gz1, x1, y1, z1);
            }

            float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
            float t22, t42;
            if (t2 < 0.0) n2 = t2 = t22 = t42 = gx2 = gy2 = gz2 = 0.0f;
            else
            {
                Perlin.gradrot(Perlin.perm[ii + i2 + Perlin.perm[jj + j2 + Perlin.perm[kk + k2]]], sin_t, cos_t, out gx2, out gy2, out gz2);
                t22 = t2 * t2;
                t42 = t22 * t22;
                n2 = t42 * Perlin.graddot(gx2, gy2, gz2, x2, y2, z2);
            }

            float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
            float t23, t43;
            if (t3 < 0.0) n3 = t3 = t23 = t43 = gx3 = gy3 = gz3 = 0.0f;
            else
            {
                Perlin.gradrot(Perlin.perm[ii + 1 + Perlin.perm[jj + 1 + Perlin.perm[kk + 1]]], sin_t, cos_t, out gx3, out gy3, out gz3);
                t23 = t3 * t3;
                t43 = t23 * t23;
                n3 = t43 * Perlin.graddot(gx3, gy3, gz3, x3, y3, z3);
            }

            /*  Add contributions from each corner to get the final noise value.
             * The result is scaled to return values in the range [- 1,1] */
            noise = 28.0f * (n0 + n1 + n2 + n3);

            /* Compute derivative, if requested by supplying non - null pointers
             * for the last three arguments */
            float dnoise_dx, dnoise_dy, dnoise_dz;

            /*  A straight, unoptimised calculation would be like:
             *     * dnoise_dx = - 8.0 * t20 * t0 * x0 * Perlin.graddot(gx0, gy0, gz0, x0, y0, z0) + t40 * gx0;
             *    * dnoise_dy = - 8.0 * t20 * t0 * y0 * Perlin.graddot(gx0, gy0, gz0, x0, y0, z0) + t40 * gy0;
             *    * dnoise_dz = - 8.0 * t20 * t0 * z0 * Perlin.graddot(gx0, gy0, gz0, x0, y0, z0) + t40 * gz0;
             *    * dnoise_dx += - 8.0 * t21 * t1 * x1 * Perlin.graddot(gx1, gy1, gz1, x1, y1, z1) + t41 * gx1;
             *    * dnoise_dy += - 8.0 * t21 * t1 * y1 * Perlin.graddot(gx1, gy1, gz1, x1, y1, z1) + t41 * gy1;
             *    * dnoise_dz += - 8.0 * t21 * t1 * z1 * Perlin.graddot(gx1, gy1, gz1, x1, y1, z1) + t41 * gz1;
             *    * dnoise_dx += - 8.0 * t22 * t2 * x2 * Perlin.graddot(gx2, gy2, gz2, x2, y2, z2) + t42 * gx2;
             *    * dnoise_dy += - 8.0 * t22 * t2 * y2 * Perlin.graddot(gx2, gy2, gz2, x2, y2, z2) + t42 * gy2;
             *    * dnoise_dz += - 8.0 * t22 * t2 * z2 * Perlin.graddot(gx2, gy2, gz2, x2, y2, z2) + t42 * gz2;
             *    * dnoise_dx += - 8.0 * t23 * t3 * x3 * Perlin.graddot(gx3, gy3, gz3, x3, y3, z3) + t43 * gx3;
             *    * dnoise_dy += - 8.0 * t23 * t3 * y3 * Perlin.graddot(gx3, gy3, gz3, x3, y3, z3) + t43 * gy3;
             *    * dnoise_dz += - 8.0 * t23 * t3 * z3 * Perlin.graddot(gx3, gy3, gz3, x3, y3, z3) + t43 * gz3;
             */
            float temp0 = t20 * t0 * Perlin.graddot(gx0, gy0, gz0, x0, y0, z0);
            dnoise_dx = temp0 * x0;
            dnoise_dy = temp0 * y0;
            dnoise_dz = temp0 * z0;
            float temp1 = t21 * t1 * Perlin.graddot(gx1, gy1, gz1, x1, y1, z1);
            dnoise_dx += temp1 * x1;
            dnoise_dy += temp1 * y1;
            dnoise_dz += temp1 * z1;
            float temp2 = t22 * t2 * Perlin.graddot(gx2, gy2, gz2, x2, y2, z2);
            dnoise_dx += temp2 * x2;
            dnoise_dy += temp2 * y2;
            dnoise_dz += temp2 * z2;
            float temp3 = t23 * t3 * Perlin.graddot(gx3, gy3, gz3, x3, y3, z3);
            dnoise_dx += temp3 * x3;
            dnoise_dy += temp3 * y3;
            dnoise_dz += temp3 * z3;
            dnoise_dx *= -8.0f;
            dnoise_dy *= -8.0f;
            dnoise_dz *= -8.0f;
            /* This corrects a bug in the original implementation */
            dnoise_dx += t40 * gx0 + t41 * gx1 + t42 * gx2 + t43 * gx3;
            dnoise_dy += t40 * gy0 + t41 * gy1 + t42 * gy2 + t43 * gy3;
            dnoise_dz += t40 * gz0 + t41 * gz1 + t42 * gz2 + t43 * gz3;
            dnoise_dx *= 28.0f; /* Scale derivative to match the noise scaling */
            dnoise_dy *= 28.0f;
            dnoise_dz *= 28.0f;

            rx = noise;
            ry = dnoise_dx;
            rz = dnoise_dy;
            rw = dnoise_dz;
        }

        public static void GenerateCurl(float x, float y, float t, out float rx, out float ry)
        {
            GenerateDerivitives(x, y, t, out float dx, out float dy, out float dz);
            rx = dz;
            ry = -dy;
        }

        private static int FastFloor(float x)
        {
            return (x > 0) ? ((int)x) : (((int)x) - 1);
        }
    }
}
