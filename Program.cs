using System;

namespace Integral{
    class Program{
        public static double Moment0(double lim1,double lim2)
            => 3 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 2 
                        - 3 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 2;
        public static double Moment1(double lim1,double lim2)
            =>  3.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 5 
                    + 9.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 4 
                    - (3.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 5 
                    + 9.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 4);
        
        public static double Moment2(double lim1,double lim2)
            =>  3.0 * Math.Pow((lim2 - 3.0 / 2),(8.0 / 3)) / 8 
                    + 9.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 5 
                    + 27.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 8 

                    - (3.0 * Math.Pow((lim1 - 3.0 / 2),(8.0 / 3)) / 8 
                    + 9.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 5 
                    + 27.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 8);

        public static double Moment3(double lim1,double lim2)
            => 3.0 * Math.Pow((lim2 - 3.0 / 2),(11.0 / 3)) / 11 
                    + 27.0 * Math.Pow((lim2 - 3.0 / 2),(8.0 / 3)) / 16 
                    + 81.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 20 
                    + 81.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 16
                        
                    - (3.0 * Math.Pow((lim1 - 3.0 / 2),(11.0 / 3)) / 11 
                    + 27.0 * Math.Pow((lim1 - 3.0 / 2),(8.0 / 3)) / 16 
                    + 81.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 20 
                    + 81.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 16);
        public static double Moment4(double lim1,double lim2)
            => 3.0 * Math.Pow((lim2 - 3.0 / 2),(14.0 / 3)) / 14 
                    + 18.0 * Math.Pow((lim2 - 3.0 / 2),(11.0 / 3)) / 11 
                    + 81.0 * Math.Pow((lim2 - 3.0 / 2),(8.0 / 3)) / 16 
                    + 81.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 10 
                    + 243.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 32

                    - (3.0 * Math.Pow((lim1 - 3.0 / 2),(14.0 / 3)) / 14 
                    + 18.0 * Math.Pow((lim1 - 3.0 / 2),(11.0 / 3)) / 11 
                    + 81.0 * Math.Pow((lim1 - 3.0 / 2),(8.0 / 3)) / 16 
                    + 81.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 10 
                    + 243.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 32);
        public static double Moment5(double lim1,double lim2)
            => 3.0 * Math.Pow((lim2 - 3.0 / 2),(17.0 / 3)) / 17 
                    + 45.0 * Math.Pow((lim2 - 3.0 / 2),(14.0 / 3)) / 28 
                    + 135.0 * Math.Pow((lim2 - 3.0 / 2),(11.0 / 3)) / 22 
                    + 405.0 * Math.Pow((lim2 - 3.0 / 2),(8.0 / 3)) / 32 
                    + 243.0 * Math.Pow((lim2 - 3.0 / 2),(5.0 / 3)) / 16 
                    + 729.0 * Math.Pow((lim2 - 3.0 / 2),(2.0 / 3)) / 64

                    - (3.0 * Math.Pow((lim1 - 3.0 / 2),(17.0 / 3)) / 17 
                    + 45.0 * Math.Pow((lim1 - 3.0 / 2),(14.0 / 3)) / 28 
                    + 135.0 * Math.Pow((lim1 - 3.0 / 2),(11.0 / 3)) / 22 
                    + 405.0 * Math.Pow((lim1 - 3.0 / 2),(8.0 / 3)) / 32 
                    + 243.0 * Math.Pow((lim1 - 3.0 / 2),(5.0 / 3)) / 16 
                    + 729.0 * Math.Pow((lim1 - 3.0 / 2),(2.0 / 3)) / 64);

        public static double[] Kardano(double[] a, double[] x){
            double cur = a[0]; 
            a[0] = a[2]; 
            a[2] = cur;
            double p = a[1] - a[0] * a[0] / 3.0; 
            double q = a[2] + 2.0 * a[0] * a[0] * a[0] / 27.0 - a[0] * a[1] / 3.0;
            double determinant = q * q / 4.0 + p * p * p / 27.0;
             
            if (determinant < 0){
                double phi = 0;
                 if (q < 0) phi = Math.Atan(2.0 * Math.Sqrt(-determinant) / (-q));
                 if (q > 0) phi = Math.Atan(2.0 * Math.Sqrt(-determinant) / (-q) + Math.PI);
                 if (q == 0) phi = Math.PI / 2.0;

                 x[0] = 2.0 * Math.Sqrt(-p / 3.0) * Math.Cos(phi / 3.0) - a[0] / 3.0;
                 x[1] = 2.0 * Math.Sqrt(-p / 3.0) * Math.Cos(phi / 3.0 + 2.0 * Math.PI / 3.0) - a[0] / 3.0;
                 x[2] = 2.0 * Math.Sqrt(-p / 3.0) * Math.Cos(phi / 3.0 + 4.0 * Math.PI / 3.0) - a[0] / 3.0;
             }
             if (determinant > 0){
                x[1] = 0;
                if ((-q) / 2.0 + Math.Pow(determinant, 1.0 / 2.0) < 0)
                    x[1] += -Math.Pow((q) / 2.0 - Math.Pow(determinant, 1.0 / 2.0), 1.0 / 3.0);
                else x[1] += Math.Pow((-q) / 2.0 + Math.Pow(determinant, 1.0 / 2.0), 1.0 / 3.0);
                if (-q / 2.0 - Math.Pow(determinant, 1.0 / 2.0) < 0)
                    x[1] += -Math.Pow(q / 2.0 + Math.Pow(determinant, 1.0 / 2.0), 1.0 / 3.0) - a[0] / 3.0;
                else x[1] += Math.Pow(-q / 2.0 - Math.Pow(determinant, 1.0 / 2.0), 1.0 / 3.0) - a[0] / 3.0;
             }
             if (determinant == 0){
                 x[0] = 2 * Math.Pow(-q / 2.0, 1.0 / 3.0) - a[0] / 3.0;
                 x[1] =  -Math.Pow(-q / 2.0, 1.0 / 3.0) - a[0] / 3.0;
             }
            return x;
        }
        static double f(double x)
            => 2 * Math.Cos(2.5 * x) * Math.Exp(x * 1.0 / 3) + 4 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + x;
        static double NewtonKotsKF(double lim1, double lim2, double step){
            double integral = 0;
            double begin = lim1;
            var k = Math.Ceiling((b - lim2) / step) + 1;
            for(int i = 1; i <= k; i++){
                double[] points;
                double[] moments;
                moments = new double[3];

                moments[0] = Moment0(lim1,lim2);
                moments[1] = Moment1(lim1,lim2);
                moments[2] = Moment2(lim1,lim2);

                points = new double[3];
                points[0] = lim1;
                points[1] = lim1 + (lim2-lim1)/2;
                points[2] = lim2;

                double[,] matrixArray = new double[3,3];
                matrixArray[0,0] = 1;
                matrixArray[0,1] = 1;
                matrixArray[0,2] = 1;
                matrixArray[1,0] = points[0];
                matrixArray[1,1] = points[1];
                matrixArray[1,2] = points[2];
                matrixArray[2,0] = points[0] * points[0];
                matrixArray[2,1] = points[1] * points[1];
                matrixArray[2,2] = points[2] * points[2];

                Matrix matrix = new Matrix(matrixArray);
                matrix.LUDecompose();

                var solvation = matrix.Solution(moments);

                for(int j = 0; j < 3; j++)
                    integral += solvation[0,j] * f(points[j]);

                lim1 = begin + (i) * step;
                lim2 = begin + (i + 1) * step;
            }
            return integral;
        }
        
        static double GaussKF(double lim1,double lim2,double step){
            double result = 0;
            double begin = lim1;
            var k = Math.Ceiling((b - lim2) / step) + 1;
            for(int i = 1; i <= k; i++){
                double[] moments = new double[6];
                moments[0] = Moment0(lim1,lim2);
                moments[1] = Moment1(lim1,lim2);
                moments[2] = Moment2(lim1,lim2);
                moments[3] = Moment3(lim1,lim2);
                moments[4] = Moment4(lim1,lim2);
                moments[5] = Moment5(lim1,lim2);

                double[] points = new double[3];
                points[0] = lim1;
                points[1] = lim1 + (lim2-lim1)/2;
                points[2] = lim2;

                double[,] matArray = new double[3,3];
                matArray[0,0] = moments[0];
                matArray[0,1] = moments[1];
                matArray[0,2] = moments[2];
                matArray[1,0] = moments[1];
                matArray[1,1] = moments[2];
                matArray[1,2] = moments[3];
                matArray[2,0] = moments[2];
                matArray[2,1] = moments[3];
                matArray[2,2] = moments[4];
                
                Matrix matrix = new Matrix(matArray);
                matrix.LUDecompose();

                double[] y = new double[3];
                y[0] = -moments[3];
                y[1] = -moments[4];
                y[2] = -moments[5];

                var solution = matrix.Solution(y);
                double[] solvation = new double[solution.GetLength(1)];
                for(int j = 0; j < solution.GetLength(1); j++)
                    solvation[j] = solution[0,j];

                var p = Kardano(solvation,points);

                Array.Sort(p);

                matrix._mat[0,0] = 1;
                matrix._mat[0,1] = 1;
                matrix._mat[0,2] = 1;
                matrix._mat[1,0] = p[0];
                matrix._mat[1,1] = p[1];
                matrix._mat[1,2] = p[2];
                matrix._mat[2,0] = p[0] * p[0];
                matrix._mat[2,1] = p[1] * p[1];
                matrix._mat[2,2] = p[2] * p[2];

                matrix.LUDecompose();
                
                var s = matrix.Solution(moments);
                double[] A = new double[s.GetLength(1)];
                for(int j = 0;j < A.Length;j++)
                    A[j] = s[0,j];

                for (int j = 0; j < 3; j++){
                    result += A[j] * f(points[j]);
                }

                lim1 = begin + (i) * step;
                lim2 = begin + (i + 1) * step;
            }
            return result;
        }
        static double a = 1.5;
        static    double b = 3.3;
        static    double alpha = 1/3;
        static    double beta = 0;
        static    double tableValue = 7.077031437995793610263911711602477164432;
        static void Main(string[] args){
            double integral = NewtonKotsKF(a,b,1);
            double metodicMistake = 0.6406;

            WriteText("Newton - Kotz :",integral,Math.Abs(tableValue - integral),metodicMistake);

            double h = Math.Ceiling((b - a) / 10);
            double begin = a;
            double step = (b - a) / 2;
            double end = 0;
            double error = 10;
            double L = 2;
            double degree = 4;
            double speed = 3;
            double integral2, integral3;
            Console.WriteLine("Сonvergence rate of a composite interpolation quadrature form :");
            while(error > 0.000001){
                end = begin + step;
                integral = NewtonKotsKF(begin,end,step);

                step = step / L;
                end = begin + step;
                integral2 = NewtonKotsKF(begin,end,step);

                step = step / L;
                end = begin + step;
                integral3 = NewtonKotsKF(begin,end,step);

                speed = -Math.Log(Math.Abs((integral3 - integral2) / (integral2 - integral))) / Math.Log(L);
                Console.WriteLine(speed);

                error = Math.Abs((integral2 - integral) / (Math.Pow(L, degree) - 1)); 
                integral += (integral2 - integral) / (1 - Math.Pow(L, -degree));
            }
            Console.WriteLine("\n");
            WriteText("Composite interpolation quadrature form :",integral,Math.Abs(tableValue - integral),error);

            double sOpt = Math.Ceiling( (b - a) / (step * L * Math.Pow((0.00001 / error), 1.0 / degree)));
            sOpt = (b - a) / sOpt;
            begin = a;
            end = a + sOpt;
            integral2 = NewtonKotsKF(begin, end, sOpt);
            error = Math.Abs((integral2 - integral) / (Math.Pow(L, degree) - 1)); 
            
            WriteText("Composite interpolation quadrature form with optimal step :",integral,Math.Abs(tableValue - integral2),error);

            begin = a;
            step = (b - a) / 10;
            error = 0.6406;

            Console.WriteLine("Convergence rate of the composite quadratic form of a Gauss :" + "\n");

            while(error > 1e-6){
                step *= L;
                end = a + step;
                integral = GaussKF(begin, end, step);

                step = step / L;
                end = a + step;
                integral2 = GaussKF(begin, end, step);

                step = step / L ;
                end = a + step;
                integral3 = GaussKF(begin, end, step);

                speed = -Math.Log(Math.Abs((integral3 - integral2) / (integral2 - integral))) / Math.Log(L);
                System.Console.Out.Write(speed + "\n");
                error = Math.Abs((integral2 - integral) / (Math.Pow(L, degree) - 1));
                integral += (integral2 - integral) / (1 - Math.Pow(L, -degree));
            }

            WriteText("Composite interpolation quadrature of a Gauss :", integral, Math.Abs(tableValue - integral), error);
            sOpt = Math.Ceiling((b - a) / (step * L * Math.Pow((0.00001 / error), 1.0 / degree)));
            sOpt = (b - a) / sOpt;
            begin = a;
            end = a + sOpt;
            integral2 = GaussKF(begin, end, sOpt);
            error = Math.Abs((integral2 - integral) / (Math.Pow(L, degree) - 1));

            WriteText("Composite interpolation quadrature of a Gauss with optimal step :", integral2, Math.Abs(tableValue - integral2), error);
        }

        static void WriteText(string s, double num){
            System.Console.Out.Write(s + "\n");
            System.Console.Out.Write(num + "\n" + "\n");
        }

        static void WriteText(string s, double num, double num2, double num3){
            System.Console.Out.Write("\n" + s + "\n");
            System.Console.Out.Write(num + "\n");
            System.Console.Out.Write("Exact error :" + "\n");
            System.Console.Out.Write(num2 + "\n");
            System.Console.Out.Write("Methodological error :" + "\n");
            System.Console.Out.Write(num3 + "\n" + "\n");
        }
    }
}
