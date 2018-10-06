using System;

namespace Noise
{ 
    public class Perlin
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

        static Perlin()
        {
            perm = new byte[permOriginal.Length];
            permOriginal.CopyTo(perm, 0);
        }

        public static int Seed
        {
            get { return seed; }
            set
            {
                if (value == 0)
                {
                    perm = new byte[permOriginal.Length];
                    permOriginal.CopyTo(perm, 0);
                }
                else
                {
                    perm = new byte[512];
                    System.Random random = new System.Random(value);
                    random.NextBytes(perm);
                }
            }
        }

        private static int seed = 0;

        /// <summary>
        /// 1D simplex noise
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static float Generate(float x)
        {
            int i0 = FastFloor(x);
            int i1 = i0 + 1;
            float x0 = x - i0;
            float x1 = x0 - 1.0f;

            float n0, n1;

            float t0 = 1.0f - x0 * x0;
            t0 *= t0;
            n0 = t0 * t0 * grad(perm[i0 & 0xff], x0);

            float t1 = 1.0f - x1 * x1;
            t1 *= t1;
            n1 = t1 * t1 * grad(perm[i1 & 0xff], x1);
            // The maximum value of this noise is 8*(3/4)^4 = 2.53125
            // A factor of 0.395 scales to fit exactly within [-1,1]
            return 0.395f * (n0 + n1);
        }

        /// <summary>
        /// 2D simplex noise
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static float Generate(float x, float y)
        {
            const float F2 = 0.366025403f; // F2 = 0.5*(sqrt(3.0)-1.0)
            const float G2 = 0.211324865f; // G2 = (3.0-Math.sqrt(3.0))/6.0

            float n0, n1, n2; // Noise contributions from the three corners

            // Skew the input space to determine which simplex cell we're in
            float s = (x + y) * F2; // Hairy factor for 2D
            float xs = x + s;
            float ys = y + s;
            int i = FastFloor(xs);
            int j = FastFloor(ys);

            float t = (float)(i + j) * G2;
            float X0 = i - t; // Unskew the cell origin back to (x,y) space
            float Y0 = j - t;
            float x0 = x - X0; // The x,y distances from the cell origin
            float y0 = y - Y0;

            // For the 2D case, the simplex shape is an equilateral triangle.
            // Determine which simplex we are in.
            int i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
            if (x0 > y0) { i1 = 1; j1 = 0; } // lower triangle, XY order: (0,0)->(1,0)->(1,1)
            else { i1 = 0; j1 = 1; }      // upper triangle, YX order: (0,0)->(0,1)->(1,1)

            // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
            // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
            // c = (3-sqrt(3))/6

            float x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
            float y1 = y0 - j1 + G2;
            float x2 = x0 - 1.0f + 2.0f * G2; // Offsets for last corner in (x,y) unskewed coords
            float y2 = y0 - 1.0f + 2.0f * G2;

            // Wrap the integer indices at 256, to avoid indexing perm[] out of bounds
            int ii = Mod(i, 256);
            int jj = Mod(j, 256);

            // Calculate the contribution from the three corners
            float t0 = 0.5f - x0 * x0 - y0 * y0;
            if (t0 < 0.0f) n0 = 0.0f;
            else
            {
                t0 *= t0;
                n0 = t0 * t0 * grad(perm[ii + perm[jj]], x0, y0);
            }

            float t1 = 0.5f - x1 * x1 - y1 * y1;
            if (t1 < 0.0f) n1 = 0.0f;
            else
            {
                t1 *= t1;
                n1 = t1 * t1 * grad(perm[ii + i1 + perm[jj + j1]], x1, y1);
            }

            float t2 = 0.5f - x2 * x2 - y2 * y2;
            if (t2 < 0.0f) n2 = 0.0f;
            else
            {
                t2 *= t2;
                n2 = t2 * t2 * grad(perm[ii + 1 + perm[jj + 1]], x2, y2);
            }

            // Add contributions from each corner to get the final noise value.
            // The result is scaled to return values in the interval [-1,1].
            return 40.0f * (n0 + n1 + n2); // TODO: The scale factor is preliminary!
        }


        public static float Generate(float x, float y, float z)
        {
            // Simple skewing factors for the 3D case
            const float F3 = 0.333333333f;
            const float G3 = 0.166666667f;

            float n0, n1, n2, n3; // Noise contributions from the four corners

            // Skew the input space to determine which simplex cell we're in
            float s = (x + y + z) * F3; // Very nice and simple skew factor for 3D
            float xs = x + s;
            float ys = y + s;
            float zs = z + s;
            int i = FastFloor(xs);
            int j = FastFloor(ys);
            int k = FastFloor(zs);

            float t = (float)(i + j + k) * G3;
            float X0 = i - t; // Unskew the cell origin back to (x,y,z) space
            float Y0 = j - t;
            float Z0 = k - t;
            float x0 = x - X0; // The x,y,z distances from the cell origin
            float y0 = y - Y0;
            float z0 = z - Z0;

            // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
            // Determine which simplex we are in.
            int i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
            int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords

            /* This code would benefit from a backport from the GLSL version! */
            if (x0 >= y0)
            {
                if (y0 >= z0)
                { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0; } // X Y Z order
                else if (x0 >= z0) { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1; } // X Z Y order
                else { i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1; } // Z X Y order
            }
            else
            { // x0<y0
                if (y0 < z0) { i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1; } // Z Y X order
                else if (x0 < z0) { i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1; } // Y Z X order
                else { i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0; } // Y X Z order
            }

            // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
            // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
            // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
            // c = 1/6.

            float x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
            float y1 = y0 - j1 + G3;
            float z1 = z0 - k1 + G3;
            float x2 = x0 - i2 + 2.0f * G3; // Offsets for third corner in (x,y,z) coords
            float y2 = y0 - j2 + 2.0f * G3;
            float z2 = z0 - k2 + 2.0f * G3;
            float x3 = x0 - 1.0f + 3.0f * G3; // Offsets for last corner in (x,y,z) coords
            float y3 = y0 - 1.0f + 3.0f * G3;
            float z3 = z0 - 1.0f + 3.0f * G3;

            // Wrap the integer indices at 256, to avoid indexing perm[] out of bounds
            int ii = Mod(i, 256);
            int jj = Mod(j, 256);
            int kk = Mod(k, 256);

            // Calculate the contribution from the four corners
            float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
            if (t0 < 0.0f) n0 = 0.0f;
            else
            {
                t0 *= t0;
                n0 = t0 * t0 * grad(perm[ii + perm[jj + perm[kk]]], x0, y0, z0);
            }

            float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
            if (t1 < 0.0f) n1 = 0.0f;
            else
            {
                t1 *= t1;
                n1 = t1 * t1 * grad(perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]], x1, y1, z1);
            }

            float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
            if (t2 < 0.0f) n2 = 0.0f;
            else
            {
                t2 *= t2;
                n2 = t2 * t2 * grad(perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]], x2, y2, z2);
            }

            float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
            if (t3 < 0.0f) n3 = 0.0f;
            else
            {
                t3 *= t3;
                n3 = t3 * t3 * grad(perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]], x3, y3, z3);
            }

            // Add contributions from each corner to get the final noise value.
            // The result is scaled to stay just inside [-1,1]
            return 32.0f * (n0 + n1 + n2 + n3); // TODO: The scale factor is preliminary!
        }

        private static byte[,] sSimplexLut = new byte[64, 4] {
            {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 0, 0, 0}, {0, 2, 3, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 2, 3, 0},
            {0, 2, 1, 3}, {0, 0, 0, 0}, {0, 3, 1, 2}, {0, 3, 2, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 3, 2, 0},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {1, 2, 0, 3}, {0, 0, 0, 0}, {1, 3, 0, 2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
            {1, 0, 2, 3}, {1, 0, 3, 2}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {2, 0, 3, 1}, {0, 0, 0, 0}, {2, 1, 3, 0},
            {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
            {2, 0, 1, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {3, 0, 1, 2}, {3, 0, 2, 1}, {0, 0, 0, 0}, {3, 1, 2, 0},
            {2, 1, 0, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {3, 1, 0, 2}, {0, 0, 0, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}
        };

        // 4D simplex noise
        public static float Generate(float x, float y, float z, float w)
        {
            const float F4 = 0.309016994f; // F4 = (Math.sqrt(5.0)- 1.0)/ 4.0
            const float G4 = 0.138196601f; // G4 = (5.0 -Math.sqrt(5.0))/ 20.0

            float n0, n1, n2, n3, n4; // Noise contributions from the five corners

            // Skew the (x,y,z,w) space to determine which cell of 24 simplices we're in
            float s = (x + y + z + w) * F4; // Factor for 4D skewing
            float xs = x + s;
            float ys = y + s;
            float zs = z + s;
            float ws = w + s;
            int i = FastFloor(xs);
            int j = FastFloor(ys);
            int k = FastFloor(zs);
            int l = FastFloor(ws);

            float t = (i + j + k + l) * G4; // Factor for 4D unskewing
            float X0 = i - t; // Unskew the cell origin back to (x,y,z,w) space
            float Y0 = j - t;
            float Z0 = k - t;
            float W0 = l - t;

            float x0 = x - X0;  // The x,y,z,w distances from the cell origin
            float y0 = y - Y0;
            float z0 = z - Z0;
            float w0 = w - W0;

            // For the 4D case, the simplex is a 4D shape I won't even try to describe.
            // To find out which of the 24 possible simplices we're in, we need to
            // determine the magnitude ordering of x0, y0, z0 and w0.
            // The method below is a good way of finding the ordering of x,y,z,w and
            // then find the correct traversal order for the simplex weíre in.
            // First, six pair - wise comparisons are performed between each possible pair
            // of the four coordinates, and the results are used to add up binary bits
            // for an integer index.
            int c1 = (x0 > y0) ? 32 : 0;
            int c2 = (x0 > z0) ? 16 : 0;
            int c3 = (y0 > z0) ? 8 : 0;
            int c4 = (x0 > w0) ? 4 : 0;
            int c5 = (y0 > w0) ? 2 : 0;
            int c6 = (z0 > w0) ? 1 : 0;
            int c = c1 + c2 + c3 + c4 + c5 + c6;

            int i1, j1, k1, l1; // The integer offsets for the second simplex corner
            int i2, j2, k2, l2; // The integer offsets for the third simplex corner
            int i3, j3, k3, l3; // The integer offsets for the fourth simplex corner

            // sSimplexLut[c] is a 4 - vector with the numbers 0, 1, 2 and 3 in some order.
            // Many values of c will never occur, since e.g. x > y > z > w makes x < z, y < w and x < w
            // impossible. Only the 24 indices which have non - zero entries make any sense.
            // We use a thresholding to set the coordinates in turn from the largest magnitude.
            // The number 3 in the "simplex" array is at the position of the largest coordinate.
            i1 = sSimplexLut[c, 0] >= 3 ? 1 : 0;
            j1 = sSimplexLut[c, 1] >= 3 ? 1 : 0;
            k1 = sSimplexLut[c, 2] >= 3 ? 1 : 0;
            l1 = sSimplexLut[c, 3] >= 3 ? 1 : 0;
            // The number 2 in the "simplex" array is at the second largest coordinate.
            i2 = sSimplexLut[c, 0] >= 2 ? 1 : 0;
            j2 = sSimplexLut[c, 1] >= 2 ? 1 : 0;
            k2 = sSimplexLut[c, 2] >= 2 ? 1 : 0;
            l2 = sSimplexLut[c, 3] >= 2 ? 1 : 0;
            // The number 1 in the "simplex" array is at the second smallest coordinate.
            i3 = sSimplexLut[c, 0] >= 1 ? 1 : 0;
            j3 = sSimplexLut[c, 1] >= 1 ? 1 : 0;
            k3 = sSimplexLut[c, 2] >= 1 ? 1 : 0;
            l3 = sSimplexLut[c, 3] >= 1 ? 1 : 0;
            // The fifth corner has all coordinate offsets = 1, so no need to look that up.

            float x1 = x0 - i1 + G4; // Offsets for second corner in (x,y,z,w) coords
            float y1 = y0 - j1 + G4;
            float z1 = z0 - k1 + G4;
            float w1 = w0 - l1 + G4;
            float x2 = x0 - i2 + 2.0f * G4; // Offsets for third corner in (x,y,z,w) coords
            float y2 = y0 - j2 + 2.0f * G4;
            float z2 = z0 - k2 + 2.0f * G4;
            float w2 = w0 - l2 + 2.0f * G4;
            float x3 = x0 - i3 + 3.0f * G4; // Offsets for fourth corner in (x,y,z,w) coords
            float y3 = y0 - j3 + 3.0f * G4;
            float z3 = z0 - k3 + 3.0f * G4;
            float w3 = w0 - l3 + 3.0f * G4;
            float x4 = x0 - 1.0f + 4.0f * G4; // Offsets for last corner in (x,y,z,w) coords
            float y4 = y0 - 1.0f + 4.0f * G4;
            float z4 = z0 - 1.0f + 4.0f * G4;
            float w4 = w0 - 1.0f + 4.0f * G4;

            // Wrap the integer indices at 256, to avoid indexing perm[] out of bounds
            int ii = i & 0xff;
            int jj = j & 0xff;
            int kk = k & 0xff;
            int ll = l & 0xff;

            // Calculate the contribution from the five corners
            float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0 - w0 * w0;
            if (t0 < 0.0f) n0 = 0.0f;
            else
            {
                t0 *= t0;
                n0 = t0 * t0 * grad(perm[ii + perm[jj + perm[kk + perm[ll]]]], x0, y0, z0, w0);
            }

            float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1 - w1 * w1;
            if (t1 < 0.0f) n1 = 0.0f;
            else
            {
                t1 *= t1;
                n1 = t1 * t1 * grad(perm[ii + i1 + perm[jj + j1 + perm[kk + k1 + perm[ll + l1]]]], x1, y1, z1, w1);
            }

            float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2 - w2 * w2;
            if (t2 < 0.0f) n2 = 0.0f;
            else
            {
                t2 *= t2;
                n2 = t2 * t2 * grad(perm[ii + i2 + perm[jj + j2 + perm[kk + k2 + perm[ll + l2]]]], x2, y2, z2, w2);
            }

            float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3 - w3 * w3;
            if (t3 < 0.0f) n3 = 0.0f;
            else
            {
                t3 *= t3;
                n3 = t3 * t3 * grad(perm[ii + i3 + perm[jj + j3 + perm[kk + k3 + perm[ll + l3]]]], x3, y3, z3, w3);
            }

            float t4 = 0.6f - x4 * x4 - y4 * y4 - z4 * z4 - w4 * w4;
            if (t4 < 0.0f) n4 = 0.0f;
            else
            {
                t4 *= t4;
                n4 = t4 * t4 * grad(perm[ii + 1 + perm[jj + 1 + perm[kk + 1 + perm[ll + 1]]]], x4, y4, z4, w4);
            }

            // Sum up and scale the result to cover the range [- 1,1]
            return 27.0f * (n0 + n1 + n2 + n3 + n4); // TODO: The scale factor is preliminary!
        }

        public static void GenerateDerivatives(float x, out float rx, out float ry)
        {
            int i0 = FastFloor(x);
            int i1 = i0 + 1;
            float x0 = x - i0;
            float x1 = x0 - 1.0f;

            float gx0 = 0.0f, gx1 = 0.0f;
            float n0, n1;
            float t20, t40, t21, t41;

            float x20 = x0 * x0;
            float t0 = 1.0f - x20;
            //  if(t0 < 0.0) t0 = 0.0; // Never happens for 1D: x0 <= 1 always
            t20 = t0 * t0;
            t40 = t20 * t20;
            gradder(perm[i0 & 0xff], out gx0);
            n0 = t40 * gx0 * x0;

            float x21 = x1 * x1;
            float t1 = 1.0f - x21;
            //  if(t1 < 0.0) t1 = 0.0; // Never happens for 1D: |x1|<= 1 always
            t21 = t1 * t1;
            t41 = t21 * t21;
            gradder(perm[i1 & 0xff], out gx1);
            n1 = t41 * gx1 * x1;

            /* Compute derivative according to:
             *  * dnoise_dx = - 8.0 * t20 * t0 * x0 * (gx0 * x0) + t40 * gx0;
             *  * dnoise_dx += - 8.0 * t21 * t1 * x1 * (gx1 * x1) + t41 * gx1;
             */
            float dnoise_dx = t20 * t0 * gx0 * x20;
            dnoise_dx += t21 * t1 * gx1 * x21;
            dnoise_dx *= -8.0f;
            dnoise_dx += t40 * gx0 + t41 * gx1;
            dnoise_dx *= 0.25f; /* Scale derivative to match the noise scaling */

            // The maximum value of this noise is 8 *(3 / 4)^4 = 2.53125
            // A factor of 0.395 would scale to fit exactly within [- 1,1], but
            // to better match classic Perlin noise, we scale it down some more.
            //#if defined SIMPLEX_DERIVATIVES_RESCALE
            //    rx = 0.3961965135 * (n0 + n1), ry = dnoise_dx;
            //#else
            rx = 0.25f * (n0 + n1);
            ry = dnoise_dx;
//#endif
        }

        public static void GenerateDerivatives(float x, float y, out float rx, out float ry, out float rz)
        {
            const float F2 = 0.366025403f;
            const float G2 = 0.211324865f;

            float n0, n1, n2; // Noise contributions from the three corners

            // Skew the input space to determine which simplex cell we're in
            float s = (x + y) * F2; // Hairy factor for 2D
            float xs = x + s;
            float ys = y + s;
            int i = FastFloor(xs);
            int j = FastFloor(ys);

            float t = (i + j) * G2;
            float X0 = i - t; // Unskew the cell origin back to (x,y) space
            float Y0 = j - t;
            float x0 = x - X0; // The x,y distances from the cell origin
            float y0 = y - Y0;

            // For the 2D case, the simplex shape is an equilateral triangle.
            // Determine which simplex we are in.
            int i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
            if (x0 > y0)
            {
                i1 = 1; j1 = 0;
            }                             // lower triangle, XY order: (0,0)->(1,0)->(1,1)
            else
            {
                i1 = 0; j1 = 1;
            }                           // upper triangle, YX order: (0,0)->(0,1)->(1,1)

            // A step of (1,0) in (i,j) means a step of (1 - c,- c) in (x,y), and
            // a step of (0,1) in (i,j) means a step of (- c,1 - c) in (x,y), where
            // c = (3 - sqrt(3))/ 6

            float x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
            float y1 = y0 - j1 + G2;
            float x2 = x0 - 1.0f + 2.0f * G2; // Offsets for last corner in (x,y) unskewed coords
            float y2 = y0 - 1.0f + 2.0f * G2;

            // Wrap the integer indices at 256, to avoid indexing perm[] out of bounds
            int ii = i & 0xff;
            int jj = j & 0xff;

            float gx0 = 0.0f, gy0 = 0.0f, gx1 = 0.0f, gy1 = 0.0f, gx2 = 0.0f, gy2 = 0.0f; /* Gradients at simplex corners */

            /* Calculate the contribution from the three corners */
            float t0 = 0.5f - x0 * x0 - y0 * y0;
            float t20, t40;
            if (t0 < 0.0f) t40 = t20 = t0 = n0 = gx0 = gy0 = 0.0f; /* No influence */
            else
            {
                gradder(perm[ii + perm[jj]], out gx0, out gy0);
                t20 = t0 * t0;
                t40 = t20 * t20;
                n0 = t40 * (gx0 * x0 + gy0 * y0);
            }

            float t1 = 0.5f - x1 * x1 - y1 * y1;
            float t21, t41;
            if (t1 < 0.0f) t21 = t41 = t1 = n1 = gx1 = gy1 = 0.0f; /* No influence */
            else
            {
                gradder(perm[ii + i1 + perm[jj + j1]], out gx1, out gy1);
                t21 = t1 * t1;
                t41 = t21 * t21;
                n1 = t41 * (gx1 * x1 + gy1 * y1);
            }

            float t2 = 0.5f - x2 * x2 - y2 * y2;
            float t22, t42;
            if (t2 < 0.0f) t42 = t22 = t2 = n2 = gx2 = gy2 = 0.0f; /* No influence */
            else
            {
                gradder(perm[ii + 1 + perm[jj + 1]], out gx2, out gy2);
                t22 = t2 * t2;
                t42 = t22 * t22;
                n2 = t42 * (gx2 * x2 + gy2 * y2);
            }

            /* Compute derivative, if requested by supplying non - null pointers
             * for the last two arguments */
            /*  A straight, unoptimised calculation would be like:
             *    * dnoise_dx = - 8.0 * t20 * t0 * x0 * (gx0 * x0 + gy0 * y0) + t40 * gx0;
             *    * dnoise_dy = - 8.0 * t20 * t0 * y0 * (gx0 * x0 + gy0 * y0) + t40 * gy0;
             *    * dnoise_dx += - 8.0 * t21 * t1 * x1 * (gx1 * x1 + gy1 * y1) + t41 * gx1;
             *    * dnoise_dy += - 8.0 * t21 * t1 * y1 * (gx1 * x1 + gy1 * y1) + t41 * gy1;
             *    * dnoise_dx += - 8.0 * t22 * t2 * x2 * (gx2 * x2 + gy2 * y2) + t42 * gx2;
             *    * dnoise_dy += - 8.0 * t22 * t2 * y2 * (gx2 * x2 + gy2 * y2) + t42 * gy2;
             */
            float temp0 = t20 * t0 * (gx0 * x0 + gy0 * y0);
            float dnoise_dx = temp0 * x0;
            float dnoise_dy = temp0 * y0;
            float temp1 = t21 * t1 * (gx1 * x1 + gy1 * y1);
            dnoise_dx += temp1 * x1;
            dnoise_dy += temp1 * y1;
            float temp2 = t22 * t2 * (gx2 * x2 + gy2 * y2);
            dnoise_dx += temp2 * x2;
            dnoise_dy += temp2 * y2;
            dnoise_dx *= -8.0f;
            dnoise_dy *= -8.0f;
            dnoise_dx += t40 * gx0 + t41 * gx1 + t42 * gx2;
            dnoise_dy += t40 * gy0 + t41 * gy1 + t42 * gy2;
            dnoise_dx *= 40.0f; /* Scale derivative to match the noise scaling */
            dnoise_dy *= 40.0f;

            // Add contributions from each corner to get the final noise value.
            // The result is scaled to return values in the interval [- 1,1].
            /*#if defined SIMPLEX_DERIVATIVES_RESCALE
                rx = 70.175438596 * (n0 + n1 + n2);     // TODO: The scale factor is preliminary!
                ry = dnoise_dx;
                rz = dnoise_dy;
            #else*/
            rx = 40.0f * (n0 + n1 + n2);
            ry = dnoise_dx;
            rz = dnoise_dy;     // TODO: The scale factor is preliminary!
                                //#endif
        }

        public static void GenerateDerivatives(float x, float y, float z, out float rx, out float ry, out float rz, out float rw)
        {
            const float F3 = 0.333333333f;
            const float G3 = 0.166666667f;

            float n0, n1, n2, n3; /* Noise contributions from the four simplex corners */
            float rnoise;          /* Return value */
            float gx0 = 0.0f, gy0 = 0.0f, gz0 = 0.0f, gx1 = 0.0f, gy1 = 0.0f, gz1 = 0.0f; /* Gradients at simplex corners */
            float gx2 = 0.0f, gy2 = 0.0f, gz2 = 0.0f, gx3 = 0.0f, gy3 = 0.0f, gz3 = 0.0f;

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

            /* Wrap the integer indices at 256, to avoid indexing perm[] out of bounds */
            int ii = i & 0xff;
            int jj = j & 0xff;
            int kk = k & 0xff;

            /* Calculate the contribution from the four corners */
            float t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
            float t20, t40;
            if (t0 < 0.0f) n0 = t0 = t20 = t40 = gx0 = gy0 = gz0 = 0.0f;
            else
            {
                gradder(perm[ii + perm[jj + perm[kk]]], out gx0, out gy0, out gz0);
                t20 = t0 * t0;
                t40 = t20 * t20;
                n0 = t40 * (gx0 * x0 + gy0 * y0 + gz0 * z0);
            }

            float t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
            float t21, t41;
            if (t1 < 0.0f) n1 = t1 = t21 = t41 = gx1 = gy1 = gz1 = 0.0f;
            else
            {
                gradder(perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]], out gx1, out gy1, out gz1);
                t21 = t1 * t1;
                t41 = t21 * t21;
                n1 = t41 * (gx1 * x1 + gy1 * y1 + gz1 * z1);
            }

            float t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
            float t22, t42;
            if (t2 < 0.0f) n2 = t2 = t22 = t42 = gx2 = gy2 = gz2 = 0.0f;
            else
            {
                gradder(perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]], out gx2, out gy2, out gz2);
                t22 = t2 * t2;
                t42 = t22 * t22;
                n2 = t42 * (gx2 * x2 + gy2 * y2 + gz2 * z2);
            }

            float t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
            float t23, t43;
            if (t3 < 0.0f) n3 = t3 = t23 = t43 = gx3 = gy3 = gz3 = 0.0f;
            else
            {
                gradder(perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]], out gx3, out gy3, out gz3);
                t23 = t3 * t3;
                t43 = t23 * t23;
                n3 = t43 * (gx3 * x3 + gy3 * y3 + gz3 * z3);
            }

            /*  Add contributions from each corner to get the final rnoise value.
             * The result is scaled to return values in the range [- 1,1] */
/*#if defined SIMPLEX_DERIVATIVES_RESCALE
    rnoise = 34.525277436f * (n0 + n1 + n2 + n3);
#else*/
            rnoise = 28.0f * (n0 + n1 + n2 + n3);
//#endif

            /* Compute derivative, if requested by supplying non - null pointers
             * for the last three arguments */
            /*  A straight, unoptimised calculation would be like:
             *     * dnoise_dx = - 8.0 * t20 * t0 * x0 * dot(gx0, gy0, gz0, x0, y0, z0) + t40 * gx0;
             *    * dnoise_dy = - 8.0 * t20 * t0 * y0 * dot(gx0, gy0, gz0, x0, y0, z0) + t40 * gy0;
             *    * dnoise_dz = - 8.0 * t20 * t0 * z0 * dot(gx0, gy0, gz0, x0, y0, z0) + t40 * gz0;
             *    * dnoise_dx += - 8.0 * t21 * t1 * x1 * dot(gx1, gy1, gz1, x1, y1, z1) + t41 * gx1;
             *    * dnoise_dy += - 8.0 * t21 * t1 * y1 * dot(gx1, gy1, gz1, x1, y1, z1) + t41 * gy1;
             *    * dnoise_dz += - 8.0 * t21 * t1 * z1 * dot(gx1, gy1, gz1, x1, y1, z1) + t41 * gz1;
             *    * dnoise_dx += - 8.0 * t22 * t2 * x2 * dot(gx2, gy2, gz2, x2, y2, z2) + t42 * gx2;
             *    * dnoise_dy += - 8.0 * t22 * t2 * y2 * dot(gx2, gy2, gz2, x2, y2, z2) + t42 * gy2;
             *    * dnoise_dz += - 8.0 * t22 * t2 * z2 * dot(gx2, gy2, gz2, x2, y2, z2) + t42 * gz2;
             *    * dnoise_dx += - 8.0 * t23 * t3 * x3 * dot(gx3, gy3, gz3, x3, y3, z3) + t43 * gx3;
             *    * dnoise_dy += - 8.0 * t23 * t3 * y3 * dot(gx3, gy3, gz3, x3, y3, z3) + t43 * gy3;
             *    * dnoise_dz += - 8.0 * t23 * t3 * z3 * dot(gx3, gy3, gz3, x3, y3, z3) + t43 * gz3;
             */
            float temp0 = t20 * t0 * (gx0 * x0 + gy0 * y0 + gz0 * z0);
            float dnoise_dx = temp0 * x0;
            float dnoise_dy = temp0 * y0;
            float dnoise_dz = temp0 * z0;
            float temp1 = t21 * t1 * (gx1 * x1 + gy1 * y1 + gz1 * z1);
            dnoise_dx += temp1 * x1;
            dnoise_dy += temp1 * y1;
            dnoise_dz += temp1 * z1;
            float temp2 = t22 * t2 * (gx2 * x2 + gy2 * y2 + gz2 * z2);
            dnoise_dx += temp2 * x2;
            dnoise_dy += temp2 * y2;
            dnoise_dz += temp2 * z2;
            float temp3 = t23 * t3 * (gx3 * x3 + gy3 * y3 + gz3 * z3);
            dnoise_dx += temp3 * x3;
            dnoise_dy += temp3 * y3;
            dnoise_dz += temp3 * z3;
            dnoise_dx *= -8.0f;
            dnoise_dy *= -8.0f;
            dnoise_dz *= -8.0f;
            dnoise_dx += t40 * gx0 + t41 * gx1 + t42 * gx2 + t43 * gx3;
            dnoise_dy += t40 * gy0 + t41 * gy1 + t42 * gy2 + t43 * gy3;
            dnoise_dz += t40 * gz0 + t41 * gz1 + t42 * gz2 + t43 * gz3;
            dnoise_dx *= 28.0f; /* Scale derivative to match the rnoise scaling */
            dnoise_dy *= 28.0f;
            dnoise_dz *= 28.0f;

            rx = rnoise;
            ry = dnoise_dx;
            rz = dnoise_dy;
            rw = dnoise_dz;
        }

        public static void GenerateCurl(float x, float y, out float rx, out float ry)
        {
            float dx, dy, dz;

            GenerateDerivatives(x, y, out dx, out dy, out dz);
            rx = dz;
            ry = -dy;
        }

        public static byte[] perm;

        private static readonly byte[] permOriginal = new byte[]
        {
            151,160,137,91,90,15,
            131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
            190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
            88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
            77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
            102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
            135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
            5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
            223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
            129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
            251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
            49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
            138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
            151,160,137,91,90,15,
            131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
            190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
            88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
            77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
            102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
            135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
            5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
            223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
            129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
            251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
            49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
            138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
        };

        private static int FastFloor(float x)
        {
            return (x > 0) ? ((int)x) : (((int)x) - 1);
        }

        private static int Mod(int x, int m)
        {
            int a = x % m;
            return a < 0 ? a + m : a;
        }

        private static float grad(int hash, float x)
        {
            int h = hash & 15;
            float grad = 1.0f + (h & 7);   // Gradient value 1.0, 2.0, ..., 8.0
            if ((h & 8) != 0) grad = -grad;         // Set a random sign for the gradient
            return (grad * x);           // Multiply the gradient with the distance
        }

        private static float grad(int hash, float x, float y)
        {
            int h = hash & 7;      // Convert low 3 bits of hash code
            float u = h < 4 ? x : y;  // into 8 simple gradient directions,
            float v = h < 4 ? y : x;  // and compute the dot product with (x,y).
            return ((h & 1) != 0 ? -u : u) + ((h & 2) != 0 ? -2.0f * v : 2.0f * v);
        }

        private static float grad(int hash, float x, float y, float z)
        {
            int h = hash & 15;     // Convert low 4 bits of hash code into 12 simple
            float u = h < 8 ? x : y; // gradient directions, and compute dot product.
            float v = h < 4 ? y : h == 12 || h == 14 ? x : z; // Fix repeats at h = 12 to 15
            return ((h & 1) != 0 ? -u : u) + ((h & 2) != 0 ? -v : v);
        }

        private static float grad(int hash, float x, float y, float z, float t)
        {
            int h = hash & 31;      // Convert low 5 bits of hash code into 32 simple
            float u = h < 24 ? x : y; // gradient directions, and compute dot product.
            float v = h < 16 ? y : z;
            float w = h < 8 ? z : t;
            return ((h & 1) != 0 ? -u : u) + ((h & 2) != 0 ? -v : v) + ((h & 4) != 0 ? -w : w);
        }

        private static void gradder(int hash, out float x)
        {
            int h = hash & 15;
            x = 1.0f + (h & 7);
            if ((h & 8) != 0) x = -x;
        }

        private static void gradder(int hash, out float gx, out float gy)
        {
            int h = hash & 7;
            gx = grad2lut[h, 0];
            gy = grad2lut[h, 1];
            return;
        }

        private static void gradder(int hash, out float gx, out float gy, out float gz)
        {
            int h = hash & 15;
            gx = grad3lut[h, 0];
            gy = grad3lut[h, 1];
            gz = grad3lut[h, 2];
            return;
        }

        private static void gradder(int hash, out float gx, out float gy, out float gz, out float gw)
        {
            int h = hash & 31;
            gx = grad4lut[h, 0];
            gy = grad4lut[h, 1];
            gz = grad4lut[h, 2];
            gw = grad4lut[h, 3];
            return;
        }

        public static void gradrot(int hash, float sin_t, float cos_t, out float gx, out float gy)
        {
            int h = hash & 7;
            float gx0 = grad2lut[h, 0];
            float gy0 = grad2lut[h, 1];
            gx = cos_t * gx0 - sin_t * gy0;
            gy = sin_t * gx0 + cos_t * gy0;
        }
        public static void gradrot(int hash, float sin_t, float cos_t, out float gx, out float gy, out float gz)
        {
            int h = hash & 15;
            float gux = grad3u[h, 0];
            float guy = grad3u[h, 1];
            float guz = grad3u[h, 2];
            float gvx = grad3v[h, 0];
            float gvy = grad3v[h, 1];
            float gvz = grad3v[h, 2];
            gx = cos_t * gux + sin_t * gvx;
            gy = cos_t * guy + sin_t * gvy;
            gz = cos_t * guz + sin_t * gvz;
        }
        public static float graddot(float gx, float gy, float x, float y)
        {
            return gx * x + gy * y;
        }
        public static float graddot(float gx, float gy, float gz, float x, float y, float z)
        {
            return gx * x + gy * y + gz * z;
        }

        private static float[,] grad2lut = new float[8, 2] {
            { -1.0f, -1.0f }, { 1.0f, 0.0f }, { -1.0f, 0.0f }, { 1.0f, 1.0f },
            { -1.0f, 1.0f }, { 0.0f, -1.0f }, { 0.0f, 1.0f }, { 1.0f, -1.0f }
        };

        private static float[,] grad3lut = new float[16, 3] {
            { 1.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 1.0f }, // 12 cube edges
            { -1.0f, 0.0f, 1.0f }, { 0.0f, -1.0f, 1.0f },
            { 1.0f, 0.0f, -1.0f }, { 0.0f, 1.0f, -1.0f },
            { -1.0f, 0.0f, -1.0f }, { 0.0f, -1.0f, -1.0f },
            { 1.0f, -1.0f, 0.0f }, { 1.0f, 1.0f, 0.0f },
            { -1.0f, 1.0f, 0.0f }, { -1.0f, -1.0f, 0.0f },
            { 1.0f, 0.0f, 1.0f }, { -1.0f, 0.0f, 1.0f }, // 4 repeats to make 16
            { 0.0f, 1.0f, -1.0f }, { 0.0f, -1.0f, -1.0f }
        };

        private static float[,] grad4lut = new float[32, 4] {
            { 0.0f, 1.0f, 1.0f, 1.0f }, { 0.0f, 1.0f, 1.0f, -1.0f }, { 0.0f, 1.0f, -1.0f, 1.0f }, { 0.0f, 1.0f, -1.0f, -1.0f }, // 32 tesseract edges
            { 0.0f, -1.0f, 1.0f, 1.0f }, { 0.0f, -1.0f, 1.0f, -1.0f }, { 0.0f, -1.0f, -1.0f, 1.0f }, { 0.0f, -1.0f, -1.0f, -1.0f },
            { 1.0f, 0.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f, -1.0f }, { 1.0f, 0.0f, -1.0f, 1.0f }, { 1.0f, 0.0f, -1.0f, -1.0f },
            { -1.0f, 0.0f, 1.0f, 1.0f }, { -1.0f, 0.0f, 1.0f, -1.0f }, { -1.0f, 0.0f, -1.0f, 1.0f }, { -1.0f, 0.0f, -1.0f, -1.0f },
            { 1.0f, 1.0f, 0.0f, 1.0f }, { 1.0f, 1.0f, 0.0f, -1.0f }, { 1.0f, -1.0f, 0.0f, 1.0f }, { 1.0f, -1.0f, 0.0f, -1.0f },
            { -1.0f, 1.0f, 0.0f, 1.0f }, { -1.0f, 1.0f, 0.0f, -1.0f }, { -1.0f, -1.0f, 0.0f, 1.0f }, { -1.0f, -1.0f, 0.0f, -1.0f },
            { 1.0f, 1.0f, 1.0f, 0.0f }, { 1.0f, 1.0f, -1.0f, 0.0f }, { 1.0f, -1.0f, 1.0f, 0.0f }, { 1.0f, -1.0f, -1.0f, 0.0f },
            { -1.0f, 1.0f, 1.0f, 0.0f }, { -1.0f, 1.0f, -1.0f, 0.0f }, { -1.0f, -1.0f, 1.0f, 0.0f }, { -1.0f, -1.0f, -1.0f, 0.0f }
        };

        private static float[,] grad3u = new float[16, 3] {
            { 1.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 1.0f }, // 12 cube edges
            { -1.0f, 0.0f, 1.0f }, { 0.0f, -1.0f, 1.0f },
            { 1.0f, 0.0f, -1.0f }, { 0.0f, 1.0f, -1.0f },
            { -1.0f, 0.0f, -1.0f }, { 0.0f, -1.0f, -1.0f },
            { 0.81649658f, 0.81649658f, 0.81649658f }, { -0.81649658f, 0.81649658f, -0.81649658f },
            { -0.81649658f, -0.81649658f, 0.81649658f }, { 0.81649658f, -0.81649658f, -0.81649658f },
            { -0.81649658f, 0.81649658f, 0.81649658f }, { 0.81649658f, -0.81649658f, 0.81649658f },
            { 0.81649658f, -0.81649658f, -0.81649658f }, { -0.81649658f, 0.81649658f, -0.81649658f }
        };

        private static float[,] grad3v = new float[16, 3] {
            { -0.81649658f, 0.81649658f, 0.81649658f }, { -0.81649658f, -0.81649658f, 0.81649658f },
            { 0.81649658f, -0.81649658f, 0.81649658f }, { 0.81649658f, 0.81649658f, 0.81649658f },
            { -0.81649658f, -0.81649658f, -0.81649658f }, { 0.81649658f, -0.81649658f, -0.81649658f },
            { 0.81649658f, 0.81649658f, -0.81649658f }, { -0.81649658f, 0.81649658f, -0.81649658f },
            { 1.0f, -1.0f, 0.0f }, { 1.0f, 1.0f, 0.0f },
            { -1.0f, 1.0f, 0.0f }, { -1.0f, -1.0f, 0.0f },
            { 1.0f, 0.0f, 1.0f }, { -1.0f, 0.0f, 1.0f }, // 4 repeats to make 16
            { 0.0f, 1.0f, -1.0f }, { 0.0f, -1.0f, -1.0f }
        };
    }
}
