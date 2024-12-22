using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Text;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
        int[] a = { 5, 2, 3, 4, 5 };
        int n = a.Length/2;
        while(n > 0)
        {
            for (int i = n;i <  a.Length;i++)
            {
                int key = a[i],j = i - n;
                while(j >= 0 && a[j] > key)
                {
                    a[j + n] = a[j];
                    j -= n;
                }
                a[j + n] = key;
            }
            n /= 2;
        }
        for (int i = 0; i < a.Length; i++)
        {
            Console.WriteLine(a[i]);
        }
    }

    #region Level 1
    public int Factorial(int n)
    {
        int a = 1;
        for (int i = 2; i <= n; i++)
        {
            a *= i;
        }
        return a;
    }
    public int Combinations(int n,int k)
    {
        return (Factorial(n) / (Factorial(k) * Factorial(n - k)));
    }
    public long Task_1_1(int n, int k)
    {
        long answer = 0;

        // code here
        if (n < k || k < 0 || n < 0)
        {
            return 0;
        }
        answer = Combinations(n, k);

        // end

        return answer;
    }
    public double GeronArea(double a,double b, double c)
    {
        double p = (a + b + c) / 2;
        return Math.Sqrt(p * (p - a) * (p - b) * (p - c));
    }
    public bool exist(double a, double b, double c)
    {
        return (a < b + c) && (b < a + c) && (c < a + b);
    }


    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;

        // code here
        double a1 = first[0], b1 = first[1], c1 = first[2];
        double a2 = second[0], b2 = second[1], c2 = second[2];
        if(!exist(a1,b1,c1) || !exist(a2, b2, c2))
        {
            return -1;
        }
        // create and use GeronArea(a, b, c);
        if (GeronArea(a1, b1, c1) > GeronArea(a2, b2, c2))
            answer = 1;
        else if (GeronArea(a2, b2, c2) > GeronArea(a1, b1, c1))
            answer = 2;
        else
            answer = 0;
        // end

        // first = 1, second = 2, equal = 0, error = -1
        return answer;
    }
    public double GetDistance(double v, double a, double t)
    {
        double S = v * t + a * t * t / 2;
        return S;
    }
    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;

        // code here
        if (time < 0 || v1 < 0 || v2 < 0)
        {
            return 0;
        }
        if (GetDistance(v1,a1,time) > GetDistance(v2, a2, time)){
            answer = 1;
        }
        else if (GetDistance(v1, a1, time) < GetDistance(v2, a2, time))
        {
            answer = 2;
        }
        else
        {
            answer = 0;
        }
        // create and use GetDistance(v, a, t); t - hours

        // end

        // first = 1, second = 2, equal = 0
        return answer;
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int answer = 0;

        // code here
        int t = 1;
        while (GetDistance(v1, a1, t) > GetDistance(v2, a2, t))
            t++;

        answer = t;
        // use GetDistance(v, a, t); t - hours

        // end

        return answer;
    }
    #endregion

    #region Level 2
    public void FindMaxIndex(int[,] matrix,out int maxI, out int maxJ)
    {
        maxI = -1;
        maxJ = -1;
        int max = int.MinValue;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i,j] > max)
                {
                    max = matrix[i,j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }
    }
    public void Task_2_1(int[,] A, int[,] B)
    {
        // code here
        int amaxI,amaxJ,bmaxI,bmaxJ;
        FindMaxIndex(A, out amaxI, out amaxJ);
        FindMaxIndex(B, out bmaxI, out bmaxJ);
        int t = A[amaxI,amaxJ];
        A[amaxI,amaxJ] = B[bmaxI,bmaxJ];
        B[bmaxI, bmaxJ] = t;
        // create and use FindMaxIndex(matrix, out row, out column);

        // end
    }

    public void Task_2_2(double[] A, double[] B)
    {
        // code here

        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!

        // end
    }
    public int FindDiagonalMaxIndex(int[,] matrix)
    {
        int index = 0;
        int max = int.MinValue;
        for(int i = 0;i < matrix.GetLength(0); i++)
        {
            if (matrix[i, i] > max)
            {
                max= matrix[i, i];
                index = i;
            }
        }
        return index;
    }
    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here
        int[,] B1 = new int[4,5];
        int[,] C1 = new int[5, 6];
        for(int i = 0; i < B1.GetLength(0); i++)
        {
            for(int j = 0; j < B1.GetLength(1); j++)
            {
                if (i < FindDiagonalMaxIndex(B)){
                    B1[i, j] = B[i, j];
                }
                else
                {
                    B1[i, j] = B[i + 1, j];
                }
            }
        }
        B = B1;
        for (int i = 0;i <  C1.GetLength(0); i++)
        {
            for (int j = 0;j < C1.GetLength(1); j++)
            {
                if (i < FindDiagonalMaxIndex(C)){
                    C1[i, j] = C[i, j];
                }
                else
                {
                    C1[i,j] = C[i + 1, j];
                }
            }
        }
        C = C1;
        //  create and use method FindDiagonalMaxIndex(matrix);

        // end
    }

    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3

        // end
    }
    public void FindMaxInColumn(int[,] matrix, out int maxI)
    {
        maxI = -1;
        int max = int.MinValue;
        for (int i = 0; i < 1; i++)
        {
            for(int j = 0; j < matrix.GetLength(0); j++)
            {
                if (matrix[j, i] > max)
                {
                    max = matrix[j, i];
                    maxI = j;
                }
            }
        }
    }
    public void Task_2_5(int[,] A, int[,] B)
    {
        // code here
        int amaxI, bmaxI;
        FindMaxInColumn(A, out amaxI);
        FindMaxInColumn(B, out bmaxI);
        for (int i = 0; i < 6; ++i)
        {
            int t = B[bmaxI, i];
            B[bmaxI, i] = A[amaxI, i];
            A[amaxI, i] = t;
        }
        // create and use FindMaxInColumn(matrix, columnIndex, out rowIndex);

        // end
    }

    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here

        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);

        // end
    }
    public int CountRowPositive(int[,] matrix, int rowIndex)
    {
        int result = 0;

        int m = matrix.GetLength(1);
        for (int j = 0; j < m; j++)
        {
            if (matrix[rowIndex, j] > 0)
                result++;
        }
        return result;
    }
    public int CountColumnPositive(int[,] matrix, int colIndex)
    {
        int result = 0;

        int n = matrix.GetLength(0);

        for (int i = 0; i < n; i++)
        {
            if (matrix[i, colIndex] > 0)
                result++;
        }
        return result;
    }
    public void Task_2_7(ref int[,] B, int[,] C)
    {
        // code here

        int maxBValue = -1;
        int maxBIndex = -1;
        for (int i = 0; i < 4; i++)
        {
            if (CountRowPositive(B, i) > maxBValue)
            {
                maxBValue = CountRowPositive(B, i);
                maxBIndex = i;
            }
        }

        int maxCValue = -1;
        int maxCIndex = -1;
        for (int j = 0; j < 6; j++)
        {
            if (CountColumnPositive(C, j) > maxCValue)
            {
                maxCValue = CountColumnPositive(C, j);
                maxCIndex = j;
            }
        }

        int[,] B1 = new int[5, 5];

        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                if (i <= maxBIndex)
                    B1[i, j] = B[i, j];
                else if (i == maxBIndex + 1)
                    B1[i, j] = C[j, maxCIndex];
                else
                    B1[i, j] = B[i - 1, j];
            }
        }

        B = B1;
        // create and use CountRowPositive(matrix, rowIndex);
        // create and use CountColumnPositive(matrix, colIndex);

        // end
    }

    public void Task_2_8(int[] A, int[] B)
    {
        // code here

        // create and use SortArrayPart(array, startIndex);

        // end
    }
    public int[] SumPositiveElementsInColumns(int[,] matrix)
    {

        int[] result = new int[matrix.GetLength(1)];

        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            int summ = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, j] > 0)
                    summ += matrix[i, j];
            }

            result[j] = summ;
        }

        return result;
    }
    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = new int[A.GetLength(1) + C.GetLength(1)];
        int[] ASumm = SumPositiveElementsInColumns(A);
        int[] CSumm = SumPositiveElementsInColumns(C);

        // code here
        for (int i = 0; i < A.GetLength(1); ++i)
        {
            answer[i] = ASumm[i];
        }
        for (int i = 0; i < C.GetLength(1); ++i)
        {
            answer[i + A.GetLength(1)] = CSumm[i];
        }
        // create and use SumPositiveElementsInColumns(matrix);

        // end

        return answer;
    }

    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        // create and use RemoveColumn(matrix, columnIndex);

        // end
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here
        int amaxI, amaxJ, bmaxI, bmaxJ;
        FindMaxIndex(A, out amaxI, out amaxJ);
        FindMaxIndex(B, out bmaxI, out bmaxJ);
        int t = A[amaxI, amaxJ];
        A[amaxI, amaxJ] = B[bmaxI, bmaxJ];
        B[bmaxI, bmaxJ] = t;
        // use FindMaxIndex(matrix, out row, out column); from Task_2_1

        // end
    }
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxColumnIndex(matrix);

        // end
    }
    public void FindMinIndex(int[,] matrix, out int minI, out int minJ)
    {
        minI = -1;
        minJ = -1;
        int min = int.MaxValue;
        for(int i = 0;i < matrix.GetLength(0); ++i)
        {
            for(int j = 0; j < matrix.GetLength(1); ++j)
            {
                if(matrix[i, j] < min)
                {
                    minI = i; minJ = j;
                    min = matrix[i, j];
                }
            }
        }
    }
    public void RemoveRow(ref int[,] matrix, int rowIndex)
    {
        int[,] a = new int[matrix.GetLength(0) - 1,matrix.GetLength(1)];
        for(int i = 0;i < a.GetLength(0); ++i)
        {
            for (int j = 0;j < a.GetLength(1); ++j)
            {
                if(i < rowIndex)
                {
                    a[i, j] = matrix[i, j];
                }
                else
                {
                    a[i, j] = matrix[i + 1, j];
                }
            }
        }
        matrix = a;
    }
    public void Task_2_13(ref int[,] matrix)
    {
        // code here
        int maxI, maxJ, minI, minJ;
        FindMaxIndex(matrix, out maxI, out minJ);
        FindMinIndex(matrix, out minI, out minJ);
        RemoveRow(ref matrix, maxI);
        if(maxI > minI)
        {
            RemoveRow(ref matrix, minI);
        }
        else if (maxI < minI)
        {
            RemoveRow(ref matrix, minI - 1);
        }
        // create and use RemoveRow(matrix, rowIndex);

        // end
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here

        // create and use SortRow(matrix, rowIndex);

        // end
    }
    public double GetAverageWithoutMinMax(int[,] matrix)
    {

        int n = matrix.GetLength(0);
        int m = matrix.GetLength(1);

        int maxI, maxJ, minI, minJ;
        FindMaxIndex(matrix, out maxI, out maxJ);
        FindMinIndex(matrix, out minI, out minJ);

        int count = 0;
        double sum = 0, average = 0;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if (!(i == maxI && j == maxJ) && !(i == minJ && j == minJ))
                {
                    count++;
                    sum += matrix[i, j];
                }
            }
        }

        if (count != 0)
            average = sum / count;

        return average;
    }
    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here
        double[] sequence = { GetAverageWithoutMinMax(A), GetAverageWithoutMinMax(B), GetAverageWithoutMinMax(C) };

        if (sequence[0] < sequence[1] && sequence[1] < sequence[2])
            answer = 1;
        else if (sequence[0] > sequence[1] && sequence[1] > sequence[2])
            answer = -1;
        else
            answer = 0;
        // create and use GetAverageWithoutMinMax(matrix);

        // end

        // 1 - increasing   0 - no sequence   -1 - decreasing
        return answer;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here

        // create and use SortNegative(array);

        // end
    }
    public int[] RowsMaxElement(int[,] matrix)
    {
        int[] a = new int[matrix.GetLength(0)];
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int mx = int.MinValue;
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (mx < matrix[i, j]) mx = matrix[i, j];
            }
            a[i] = mx;
        }
        return a;

    }
    public void SortRowsByMaxElement(ref int[,] A, ref int[] a)
    {
        int n = A.GetLength(0), m = A.GetLength(1);
        for (int i = 0; i < n - 1; i++)
        {
            for (int j = 0; j < n - i - 1; j++)
            {
                if (a[j] < a[j + 1])
                {
                    int buf = a[j + 1];
                    a[j + 1] = a[j];
                    a[j] = buf;
                    for (int k = 0; k < m; k++)
                    {
                        int temp = A[j + 1, k];
                        A[j + 1, k] = A[j, k];
                        A[j, k] = temp;
                    }
                }
            }
        }
    }
    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here
        int[] a = RowsMaxElement(A), b = RowsMaxElement(B);
        SortRowsByMaxElement(ref A, ref a);
        SortRowsByMaxElement(ref B, ref b);
        // create and use SortRowsByMaxElement(matrix);

        // end
    }

    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        // create and use SortDiagonal(matrix);

        // end
    }

    public void Task_2_19(ref int[,] matrix)
    {
        // code here
        while (true)
        {
            bool fl = false;
            int n = matrix.GetLength(0);
            int m = matrix.GetLength(1);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    if (matrix[i, j] == 0)
                    {
                        RemoveRow(ref matrix, i);
                        fl = true;
                        break;
                    }
                }
                if (fl) break;
            }
            if (!fl) break;
        }
        // use RemoveRow(matrix, rowIndex); from 2_13

        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }
    public int[] CreateArrayFromMins(int[,] matrix)
    {
        int n = matrix.GetLength(0);
        int m = matrix.GetLength(1);
        int[] array = new int[n];
        for (int i = 0; i < n; i++)
        {
            int min = int.MaxValue;
            for (int j = i; j < m; j++)
            {
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                }
            }
            array[i] = min;
        }
        return array;
    }
        public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;

        // code here
        answerA = CreateArrayFromMins(A);
        answerB = CreateArrayFromMins(B);
        // create and use CreateArrayFromMins(matrix);

        // end
    }

    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here

        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);

        // end
    }
    public void Matrix_to_array(double[,] A, ref double[] arr)
    {
        int n = A.GetLength(0), m = A.GetLength(1), k = 0;
        arr = new double[n * m];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                arr[k++] = A[i, j];
            }
        }
    }

    public void rock_sort(ref double[] a)
    {
        int n = a.Length;
        for (int i = 0; i < n - 1; i++)
        {
            for (int j = 0; j < n - i - 1; j++)
            {
                if (a[j] < a[j + 1])
                {
                    double t = a[j + 1];
                    a[j + 1] = a[j];
                    a[j] = t;
                }
            }
        }
    }

    public void MatrixValuesChange(ref double[,] A, double[] mxA)
    {
        int n = A.GetLength(0), m = A.GetLength(1);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if (mxA[0] == A[i, j] || mxA[1] == A[i, j] || mxA[2] == A[i, j] || mxA[3] == A[i, j] || mxA[4] == A[i, j])
                {
                    if (A[i, j] > 0)
                    {
                        A[i, j] *= 2; 
                    }
                    else
                    {
                        A[i, j] *= 0.5; 
                    }
                }
                else
                {
                    if (A[i, j] > 0)
                    {
                        A[i, j] *= 0.5;
                    }
                    else
                    {
                        A[i, j] *= 2;
                    }
                }
            }
        }
    }
    public void Task_2_23(double[,] A, double[,] B)
    {
        // code here
        double[] mxAt = null, mxBt = null;

        Matrix_to_array(A, ref mxAt);
        Matrix_to_array(B, ref mxBt);

        rock_sort(ref mxAt);
        rock_sort(ref mxBt);

        double[] mxA = mxAt.Length >= 5 ? new double[5] : new double[mxAt.Length], mxB = mxBt.Length >= 5 ? new double[5] : new double[mxBt.Length];
        
        for (int i = 0; i < mxA.Length; i++)
        {
            mxA[i] = mxAt[i];
        }
        for (int i = 0; i < mxB.Length; i++)
        {
            mxB[i] = mxBt[i];
        }

        MatrixValuesChange(ref A, mxA);
        MatrixValuesChange(ref B, mxB);
        // create and use MatrixValuesChange(matrix);

        // end
    }

    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);

        // end
    }
    public int CountRowNegative(int[,] matrix, int rowIndex)
    {
        int result = 0;

        int m = matrix.GetLength(1);
        for (int j = 0; j < m; j++)
        {
            if (matrix[rowIndex, j] < 0)
                result++;
        }
        return result;
    }

    public int FindMaxNegativeRow(int[,] matrix)
    {
        int n = matrix.GetLength(0);

        int maxValue = -1;
        int maxRowIndex = -1;

        for (int i = 0; i < n; i++)
        {
            if (CountRowNegative(matrix, i) > maxValue)
            {
                maxValue = CountRowNegative(matrix, i);
                maxRowIndex = i;
            }
        }

        return maxRowIndex;
    }
    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;

        // code here
        maxA = FindMaxNegativeRow(A);
        maxB = FindMaxNegativeRow(B);
        // create and use FindRowWithMaxNegativeCount(matrix);
        // in FindRowWithMaxNegativeCount create and use CountNegativeInRow(matrix, rowIndex); like in 2_22

        // end
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }
    public int FindRowMaxIndex(int[,] matrix, int rowIndex)
    {
        int m = matrix.GetLength(1);
        int max = int.MinValue;
        int index = -1;
        for(int i = 0; i < m; i++)
        {
            if (matrix[rowIndex, i] > max)
            {
                max = matrix[rowIndex, i];
                index = i;
            }
        }
        return index;
    }
    public void ReplaceMaxElementOdd(int[,] matrix)
    {
        for (int i = 0;i < matrix.GetLength(0); i+=2)
        {
            int maxindex = FindRowMaxIndex(matrix, i);
            for (int j = 0;j < matrix.GetLength(1); j++)
            {
                if (matrix[i,j] == matrix[i,maxindex])
                {
                    matrix[i, j] *= j + 1;
                }
            }
        }
    }
    public void ReplaceMaxElementEven(int[,] matrix)
    {
        for (int i = 1;i < matrix.GetLength(0); i += 2)
        {
            int maxindex = FindRowMaxIndex(matrix, i);
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if(matrix[i, j] == matrix[i, maxindex])
                {
                    matrix[i, j] = 0;
                }
            }
        }
    }
    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here
        ReplaceMaxElementEven(A);
        ReplaceMaxElementEven(B);
        ReplaceMaxElementOdd(A);
        ReplaceMaxElementOdd(B);
        // create and use FindRowMaxIndex(matrix, rowIndex, columnIndex);
        // create and use ReplaceMaxElementOdd(matrix, row, column);
        // create and use ReplaceMaxElementEven(matrix, row, column);

        // end
    }

    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search

        // end
    }
    #endregion

    #region Level 3
    public delegate double SumFunction(int i, double x, ref int change);
    public delegate double YFunction(double x);
    public double s1Term(int i, double x, ref int iFactorial)
    {
        if (i > 0)
            iFactorial *= i;

        return Math.Cos(i * x) / iFactorial;
    }
    public double s2Term(int i, double x, ref int sign)
    {
        sign *= -1;
        return sign * Math.Cos(i * x) / (i * i);
    }
    public double y3_1_1(double x)
    {
        return Math.Exp(Math.Cos(x)) * Math.Cos(Math.Sin(x));
    }
    public double y3_1_2(double x)
    {
        return ((x * x) - Math.PI * Math.PI / 3) / 4;
    }
    public double CalculateSum(SumFunction sumFunction, double x, int i)
    {
        double epsilon = 0.0001, sum = 0;
        int change = 1;
        double curTerm = sumFunction(i, x, ref change);

        while (Math.Abs(curTerm) > epsilon)
        {
            sum += curTerm;
            curTerm = sumFunction(++i, x, ref change);
        }
        return sum;
    }
    public void GetSumAndY(SumFunction sFunction, YFunction yFunction, double a, double b, double h, double[,] SumAndY, int startI = 0)
    {
        for (int i = 0; i < (b - a) / h + 1; i++)
        {
            double x = a + i * h;

            double sum = CalculateSum(sFunction, x, startI);
            double y = yFunction(x);

            SumAndY[i, 0] = sum;
            SumAndY[i, 1] = y;
        }
    }
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        // code here
        double a1 = 0.1, b1 = 1, h1 = 0.1;
        firstSumAndY = new double[(int)((b1 - a1) / h1) + 1, 2];
        GetSumAndY(s1Term, y3_1_1, a1, b1, h1, firstSumAndY);


        double a2 = Math.PI / 5, b2 = Math.PI, h2 = Math.PI / 25;
        secondSumAndY = new double[(int)((b2 - a2) / h2) + 1, 2];
        GetSumAndY(s2Term, y3_1_2, a2, b2, h2, secondSumAndY, 1);
        // create and use public delegate SumFunction(x) and public delegate YFunction(x)
        // create and use method GetSumAndY(sFunction, yFunction, a, b, h);
        // create and use 2 methods for both functions calculating at specific x

        // end
    }

    public void Task_3_2(int[,] matrix)
    {
        // SortRowStyle sortStyle = default(SortRowStyle); - uncomment me

        // code here

        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting

        // end
    }
    public delegate void SwapDirection(double[] array);

    public void SwapRight(double[] array)
    {
        for (int i = 0; i < array.Length - 1; i += 2)
        {
            (array[i], array[i + 1]) = (array[i + 1], array[i]);
        }
    }

    public void SwapLeft(double[] array)
    {
        for (int i = array.Length - 1; i > 0; i -= 2)
        {
            (array[i], array[i - 1]) = (array[i - 1], array[i]);
        }
    }

    public double GetSum(double[] array)
    {
        double sum = 0;

        for (int i = 1; i < array.Length; i += 2)
        {
            sum += array[i];
        }

        return sum;
    }
    public double Task_3_3(double[] array)
    {
        double answer = 0;
         SwapDirection swapper = default(SwapDirection);

        // code here
        double sum = 0;

        foreach (double num in array)
            sum += num;

        double average = sum / array.Length;

        swapper = (array[0] > average) ? SwapRight : SwapLeft;

        swapper(array);

        answer = GetSum(array);
        // create and use public delegate SwapDirection(array);
        // create and use methods SwapRight(array) and SwapLeft(array)
        // create and use method GetSum(array, start, h) that sum elements with a negative indexes
        // change method in variable swapper in the if/else and than use swapper(matrix)

        // end

        return answer;
    }

    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here

        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end

        return answer;
    }
    public int CountSignFlips(YFunction yFunction, double a, double b, double h)
    {

        int count = 0;

        for (double x = a + h; x <= b; x += h)
        {
            double prev = yFunction(x - h), cur = yFunction(x);
            if ((prev >= 0 && cur < 0) || (prev <= 0 && cur > 0) || (x == b && cur == 0)) 
                count++;
        }

        return count;
    }

    public double y3_5_1(double x)
    {
        return x * x - Math.Sin(x);
    }

    public double y3_5_2(double x)
    {
        return Math.Exp(x) - 1;
    }
    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;

        // code here
        double a1 = 0, b1 = 2, h1 = 0.1;
        func1 = CountSignFlips(y3_5_1, a1, b1, h1);

        double a2 = -1, b2 = 1, h2 = 0.2;
        func2 = CountSignFlips(y3_5_2, a2, b2, h2);
        // use public delegate YFunction(x, a, b, h) from Task_3_1
        // create and use method CountSignFlips(YFunction, a, b, h);
        // create and use 2 methods for both functions

        // end
    }

    public void Task_3_6(int[,] matrix)
    {
        // code here

        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);

        // end
    }
    public delegate int CountPositive(int[,] matrix, int index);
    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here
        CountPositive countPositive = CountRowPositive;
        int maxB = -1;
        int maxBIndex = -1;
        for (int i = 0; i < 4; i++)
        {
            if (countPositive(B, i) > maxB)
            {
                maxB = countPositive(B, i);
                maxBIndex = i;
            }
        }
        countPositive = CountColumnPositive;
        int maxC = -1;
        int maxCIndex = -1;
        for (int j = 0; j < 6; j++)
        {
            if (countPositive(C, j) > maxC)
            {
                maxC = countPositive(C, j);
                maxCIndex = j;
            }
        }

        int[,] B1 = new int[5, 5];

        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                if (i <= maxBIndex)
                    B1[i, j] = B[i, j];
                else if (i == maxBIndex + 1)
                    B1[i, j] = C[j, maxCIndex];
                else
                    B1[i, j] = B[i - 1, j];
            }
        }

        B = B1;
        // create and use public delegate CountPositive(matrix, index);
        // use CountRowPositive(matrix, rowIndex) from Task_2_7
        // use CountColumnPositive(matrix, colIndex) from Task_2_7
        // create and use method InsertColumn(matrixB, CountRow, matrixC, CountColumn);

        // end
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // FindIndex searchArea = default(FindIndex); - uncomment me

        // code here

        // create and use public delegate FindIndex(matrix);
        // create and use method FindMaxBelowDiagonalIndex(matrix);
        // create and use method FindMinAboveDiagonalIndex(matrix);
        // use RemoveColumn(matrix, columnIndex) from Task_2_10
        // create and use method RemoveColumns(matrix, findMaxBelowDiagonalIndex, findMinAboveDiagonalIndex)

        // end
    }
    public delegate void FindElementDelegate(int[,] matrix, out int foundI, out int foundJ);
    public void RemoveRows(ref int[,] matrix, FindElementDelegate findElementDelegate1, FindElementDelegate findElementDelegate2)
    {
        int i1, j1, i2, j2;
        findElementDelegate1(matrix, out i1, out j1);
        findElementDelegate2(matrix, out i2, out j2);

        RemoveRow(ref matrix, i1);

        if (i2 < i1)
            RemoveRow(ref matrix, i2);
        else if (i2 > i1)
            RemoveRow(ref matrix, i2 - 1);
    }
    public void Task_3_13(ref int[,] matrix)
    {
        // code here
        RemoveRows(ref matrix, FindMaxIndex, FindMinIndex);
        // use public delegate FindElementDelegate(matrix) from Task_3_6
        // create and use method RemoveRows(matrix, findMax, findMin)

        // end
    }

    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);

        // end
    }
    public delegate void ReplaceMaxElement(int[,] matrix);

    public void EvenOddRowsTransform(int[,] matrix, ReplaceMaxElement replaceMaxElement1, ReplaceMaxElement replaceMaxElement2)
    {
        replaceMaxElement1(matrix);
        replaceMaxElement2(matrix);
    }
    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here
        EvenOddRowsTransform(A, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        EvenOddRowsTransform(B, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        // create and use public delegate ReplaceMaxElement(matrix, rowIndex, max);
        // use ReplaceMaxElementOdd(matrix) from Task_2_27
        // use ReplaceMaxElementEven(matrix) from Task_2_27
        // create and use method EvenOddRowsTransform(matrix, replaceMaxElementOdd, replaceMaxElementEven);

        // end
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);

        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);

        // end
    }
    #endregion
    #region bonus part
    public delegate void MatrixConverter(double[,] matrix);

    public void ToUpperTriangular(double[,] matrix)
    {
        int n = matrix.GetLength(0);
        int m = matrix.GetLength(1);

        for (int j = 0; j < m; j++)
        {
            for (int i = j + 1; i < n; i++)
            {
                double coef = -(matrix[i, j] / matrix[j, j]);

                matrix[i, j] = 0;

                for (int k = j + 1; k < m; k++)
                {
                    matrix[i, k] += matrix[j, k] * coef;
                }
            }
        }
    }

    public void ToLowerTriangular(double[,] matrix)
    {
        int n = matrix.GetLength(0);
        int m = matrix.GetLength(1);

        for (int j = m - 1; j >= 0; j--)
        {
            for (int i = j - 1; i >= 0; i--)
            {
                double coef = -(matrix[i, j] / matrix[j, j]);

                matrix[i, j] = 0;

                for (int k = j - 1; k >= 0; k--)
                {
                    matrix[i, k] += matrix[j, k] * coef;
                }
            }
        }
    }

    public void ToLeftDiagonal(double[,] matrix)
    {
        ToUpperTriangular(matrix);
        ToLowerTriangular(matrix);
    }

    public void ToRightDiagonal(double[,] matrix)
    {
        ToLowerTriangular(matrix);
        ToUpperTriangular(matrix);
    }
    public double[,] Task_4(double[,] matrix, int index)
    {
        // MatrixConverter[] mc = new MatrixConverter[4]; - uncomment me

        // code here
        MatrixConverter[] mc = { ToUpperTriangular, ToLowerTriangular, ToLeftDiagonal, ToRightDiagonal };

        mc[index](matrix);
        // create public delegate MatrixConverter(matrix);
        // create and use method ToUpperTriangular(matrix);
        // create and use method ToLowerTriangular(matrix);
        // create and use method ToLeftDiagonal(matrix); - start from the left top angle
        // create and use method ToRightDiagonal(matrix); - start from the right bottom angle

        // end

        return matrix;
    }
    #endregion
}
