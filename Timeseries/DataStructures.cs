using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DataStructures
{
    [Serializable]
    public class Series
    {
        public Series()
        {
            dateTimes = new List<DateTime>();
            yValues = new List<double>();
            normalized_yValues = new List<double>();
            
            frequencies = new List<double>();
            fftValues = new List<double>();
        }

        public Series(List<DateTime> dateTimes, List<double> yValues)
        {
            this.dateTimes = dateTimes;
            this.yValues = yValues;
            this.timeLength = yValues.Count;

            bool alreadyNormalized = (Math.Abs(mean) < 1.0e-3) && (Math.Abs(standardDeviation - 1) < 1.0e-3);
            if (!alreadyNormalized)
            {
                this.normalized_yValues = yValues;
            }
            else
            {
                this.normalized_yValues = yValues;
            }
            CalculatePSD();
            this.frequencyLength = fftValues.Count;
        }

        public void Normalize()
        {
            int N = 0;
            List<double> data = null;
            List<double> normalized_data = null;
            N = this.yValues.Count;
            normalized_data = new List<double>(normalized_yValues);
            data = new List<double>(yValues);

            // calculate mean and standard deviation for the timeseries                
            mean = 0.0;
            double sequenceValuesSum = 0.0;

            for (int i = 0; i < N; i++)
            {
                sequenceValuesSum = sequenceValuesSum + data[i];
            }
            mean = sequenceValuesSum / N;

            standardDeviation = 0.0;
            double sequenceDevSquare = 0.0;
            for (int i = 0; i < N; i++)
            {
                sequenceDevSquare = sequenceDevSquare + Math.Pow((data[i] - mean), 2);
            }
            standardDeviation = Math.Sqrt(sequenceDevSquare / N);

            // alreadyNormalized the timeseries if stdev is not zero (if timeseries fValues are not all the same)
            for (int i = 0; i < N; i++)
            {
                double normalizedValue = 0.0;
                if (standardDeviation != 0.0)
                {
                    normalizedValue = (data[i] - mean) / standardDeviation;
                }
                else
                {
                    normalizedValue = data[i];
                }

                normalized_data.Add(normalizedValue);
            }

            normalized_yValues = new List<double>(normalized_data);
        }

        public Timeseries toTimeseries()
        {
            return new Timeseries(this.dateTimes, this.yValues, this.normalized_yValues, this.mean, this.standardDeviation);
        }

        public FrequencySeries toFrequencySeries()
        {
            return new FrequencySeries(this.frequencies, this.yValues);
        }

        public void CalculatePSD()
        {
            this.frequencies.Clear();

            int n = this.timeLength;
            double[] xArray = this.normalized_yValues.ToArray<double>();
            alglib.complex[] fft;

            alglib.fftr1d(xArray, out fft);

            double[] PSD = new double[n / 2];
            for (int i = 0; i < fft.Length / 2; i++)
            {
                frequencies.Add(i);
                PSD[i] = Math.Sqrt((fft[i].x * fft[i].x) + (fft[i].y * fft[i].y));
            }

            fftValues = PSD.ToList();
        }

        public List<DateTime> dateTimes = new List<DateTime>();
        public List<double> yValues = new List<double>();
        public List<double> normalized_yValues = new List<double>();
        public List<double> frequencies = new List<double>();
        public List<double> fftValues = new List<double>();
        public int timeLength;
        public int frequencyLength;
        public double mean;
        public double standardDeviation;
    }

    [Serializable]
    public class Timeseries
    {
        public Timeseries()
        {
            dateTimes = new List<DateTime>();
            yValues = new List<double>();
            normalized_yValues = new List<double>();
        }

        public Timeseries(List<DateTime> dateTimes, List<double> yValues, List<double> normalized_yValues, double mean, double standardDeviation)
        {
            this.dateTimes = dateTimes;
            this.yValues = yValues;
            this.normalized_yValues = normalized_yValues;
            this.mean = mean;
            this.standardDeviation = standardDeviation;

            this.length = yValues.Count;
        }

        public List<DateTime> dateTimes = new List<DateTime>();
        public List<double> yValues = new List<double>();
        public List<double> normalized_yValues = new List<double>();

        public int length;
        public double mean;
        public double standardDeviation;
    }

    [Serializable]
    public class FrequencySeries
    {
        public FrequencySeries()
        {
            xValues = new List<double>();
            fValues = new List<double>();
        }

        public FrequencySeries(List<double> xValues, List<double> values)
        {
            this.xValues = xValues;
            this.fValues = values;

            this.length = values.Count;
        }

        public List<double> xValues;
        public List<double> fValues;
        public int length;
    }

    [Serializable]
    public class Location
    {
        public Location(string name)
        {
            this.name = name;
        }

        public string name;
        public List<Timeseries> timeseries = new List<Timeseries>();
        public List<int> timeWindowClusterNumbers = new List<int>();
        public List<bool> isAnomalyDetected = new List<bool>();
    }

    [Serializable]
    public class Pair
    {
        public Pair(int v1, int v2)
        {
            value1 = v1;
            value2 = v2;
        }

        public int value1;
        public int value2;
    }
}
