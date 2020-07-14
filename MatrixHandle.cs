using System;

namespace Extensions
{
    static class MatrixHandle
    {
        public static double[,] GetTransparent(this double[,] matrix)
        {
            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            var newMatrix = new double[nCols, nRows];

            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    newMatrix[j, i] = matrix[i, j];
                }
            }

            return newMatrix;
        }

        public static double[,] GetReverce(this double[,] matrix)
        {
            double determinant = matrix.GetDeterminant();
            if (Math.Abs(determinant) < 1.0e-6)
                return null;

            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);
            var newMatrix = new double[nRows, nCols];
            matrix.ProcessFunctionOverData((i, j) =>
            {
                newMatrix[i, j] = ((i + j) % 2 == 1 ? -1 : 1) *
                (matrix.CreateMatrixWithoutColumn(j).CreateMatrixWithoutRow(i).GetDeterminant()) / determinant;
            });

            return newMatrix.GetTransparent();
        }

        public static double GetDeterminant(this double[,] matrix)
        {
            double determine = 0.0d;

            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            if (nRows != nCols)
                return -1;

            determine = CalculateDeterminator(matrix);

            return determine;
        }

        public static double CalculateDeterminator(this double[,] matrix)
        {
            int nRows = matrix.GetLength(0);
            double determine = 0.0;

            if (nRows <= 0)
            {
                return -1;
            }
            else if (nRows == 1)
            {
                return matrix[0, 0];
            }
            else if (nRows == 2)
            {
                return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
            }
            else
            {
                for (var j = 0; j < nRows; j++)
                {
                    determine += (j % 2 == 1 ? 1 : -1) * matrix[1, j] *
                        matrix.CreateMatrixWithoutColumn(j).
                        CreateMatrixWithoutRow(1).CalculateDeterminator();
                }
                return determine;
            }
        }

        public static double[,] CreateMatrixWithoutColumn(this double[,] matrix, int column)
        {
            if (column < 0 || column >= matrix.GetLength(1))
            {
                return matrix;
            }

            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            var result = new double[nRows, nCols - 1];
            result.ProcessFunctionOverData((i, j) =>
                result[i, j] = j < column ? matrix[i, j] : matrix[i, j + 1]);
            return result;
        }

        public static double[,] CreateMatrixWithoutRow(this double[,] matrix, int row)
        {
            if (row < 0 || row >= matrix.GetLength(0))
            {
                return matrix;
            }

            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            var result = new double[nRows - 1, nCols];
            result.ProcessFunctionOverData((i, j) =>
                result[i, j] = i < row ? matrix[i, j] : matrix[i + 1, j]);
            return result;
        }

        public static void ProcessFunctionOverData(this double[,] matrix, Action<int, int> func)
        {
            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            for (var i = 0; i < nRows; i++)
            {
                for (var j = 0; j < nCols; j++)
                {
                    func(i, j);
                }
            }
        }

        public static double[,] MultiplyBy(this double[,] matrix, double value)
        {
            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    matrix[i, j] *= value;
                }
            }

            return matrix;
        }

        public static double[,] MultiplyBy(this double[,] matrix, double[,] multiplyMatrix)
        {
            if (matrix.GetLength(1) != multiplyMatrix.GetLength(0))
                return null;

            var newMatrix = new double[matrix.GetLength(0), multiplyMatrix.GetLength(1)];

            newMatrix.ProcessFunctionOverData((i, j) =>
            {
                for (var k = 0; k < matrix.GetLength(1); k++)
                {
                    newMatrix[i, j] += matrix[i, k] * multiplyMatrix[k, j];
                }
            });

            return newMatrix;
        }

        public static double[,] MultiplyBy(this double[,] matrix, double[] multiplyMatrix)
        {
            if (matrix.GetLength(1) != multiplyMatrix.GetLength(0))
                return null;

            var newMatrix = new double[matrix.GetLength(0), multiplyMatrix.GetLength(1)];

            newMatrix.ProcessFunctionOverData((i, j) =>
            {
                for (var k = 0; k < matrix.GetLength(1); k++)
                {
                    newMatrix[i, j] += matrix[i, k] * multiplyMatrix[k];
                }
            });

            return newMatrix;
        }

        public static double[,] NormalizeCoefficients(this double[,] matrix)
        {
            double normal = 0.0;

            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    normal += matrix[i, j] * matrix[i, j];
                }
            }

            return matrix.MultiplyBy(Math.Sqrt(normal));
        }

        public static double[,] SummWith(this double[,] matrix, double[,] otherMatrix)
        {
            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            if (nRows != otherMatrix.GetLength(0) || nCols != otherMatrix.GetLength(1))
                return null;
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    matrix[i, j] += otherMatrix[i, j];
                }
            }

            return matrix;
        }

        public static double GetMaxValue(this double[,] matrix)
        {
            double max = double.MinValue;

            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    if (matrix[i, j] > max)
                    {
                        max = matrix[i, j];
                    }
                }
            }

            return max;
        }

        public static double GetMinValue(this double[,] matrix)
        {
            double min = double.MaxValue;

            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    if (matrix[i, j] < min)
                    {
                        min = matrix[i, j];
                    }
                }
            }

            return min;
        }

        public static bool CheckEqual(this double[,] matrix, double[,] otherMatrix)
        {
            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            if (nRows != otherMatrix.GetLength(0) || nCols != otherMatrix.GetLength(1))
                return false;

            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    if (matrix[i, j] != otherMatrix[i, j])
                        return false;
                }
            }

            return true;
        }

        public static double[] ToArray(this double[,] matrix)
        {
            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);

            double[] array = new double[nRows * nCols];

            int k = 0;
            
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    array[k] = matrix[i, j];
                    ++k;
                }
            }

            return array;
        }

        public static bool CheckEqual(this double[] array, double[] other)
        {
            if (array.Length != other.Length)
                return false;

            for (int i = 0; i < array.Length; i++)
            {
                if (array[i] != other[i])
                    return false;
            }

            return true;
        }

        public static double[,] SummToAll(this double[,] matrix, double value)
        {
            int nRows = matrix.GetLength(0);
            int nCols = matrix.GetLength(1);
            
            for (int i = 0; i < nRows; i++)
            {
                for (int j = 0; j < nCols; j++)
                {
                    matrix[i, j] += value;
                }
            }

            return matrix;
        }
    }
}
