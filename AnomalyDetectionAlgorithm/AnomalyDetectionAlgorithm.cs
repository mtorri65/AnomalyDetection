using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using DataManagement;
using DataStructures;
using ConfigurationParameters;

namespace AnomalyDetection
{
    [Serializable]
    public class AnomalyDetectionAlgorithm
    {
        public AnomalyDetectionAlgorithm(List<string> locationNames, int timeWindowNumber)
        {
            this.timeWindowNumber = timeWindowNumber;
            
            foreach(string locationName in locationNames)
            {
                Location location = new Location(locationName);
                locations.Add(location);
            }
        }

        public void Detect(Parameters parameters, string locationName)
        {
            Location location = locations.Find(x => x.name == locationName);

            Timeseries thisTimeWindowSample = DataManager.ReadTimeSeries(location.name, timeWindowNumber);
            if (timeWindowNumber == 1)
            {
                List<List<double>> DTWtimeWindowDistances = new List<List<double>>();
                DTWtimeWindowDistances.Add(new List<double>(new double[1]));
                DTWdistancesEvolution.Add(location.name, DTWtimeWindowDistances);

                List<Timeseries> lts = new List<Timeseries>();
                lts.Add(thisTimeWindowSample);
                timeSamples.Add(location.name, lts);

                return;
            }

            // add one extra column (and one extra row) to the Distance Matrix by including this timeWindowNumber's timeseries in it
            for (int row = 0; row < DTWdistancesEvolution[location.name].Count; row++)
            {
                DTWdistancesEvolution[location.name][row].Add(0.0);
            }
            DTWdistancesEvolution[location.name].Add((new double[DTWdistancesEvolution[location.name].Count + 1]).ToList());

            for (int i = 0; i < DTWdistancesEvolution[location.name].Count; i++)
            {
                if (i < DTWdistancesEvolution[location.name].Count - 1)
                {
                    Timeseries sampleShifted = ShiftTimeseries(timeSamples[location.name][i], thisTimeWindowSample);
                    List<DateTime> aggregateDatetimes = CreateDatetimesAggregate(timeSamples[location.name][i], sampleShifted);
                    Timeseries sampleStretched1 = StretchTimeseries(timeSamples[location.name][i], aggregateDatetimes);
                    Timeseries sampleStretched2 = StretchTimeseries(sampleShifted, aggregateDatetimes);

                    DTWdistancesEvolution[location.name][i][DTWdistancesEvolution[location.name].Count - 1] = CalculateDTWDistance(sampleStretched1, sampleStretched2);
                    DTWdistancesEvolution[location.name][DTWdistancesEvolution[location.name].Count - 1][i] = DTWdistancesEvolution[location.name][i][DTWdistancesEvolution[location.name].Count - 1];
                }
                else
                {
                    DTWdistancesEvolution[location.name][i][DTWdistancesEvolution[location.name].Count - 1] = 0.0; // diagonal term
                }
            }
            timeSamples[location.name].Add(thisTimeWindowSample);
            CreateClusters(location, parameters);
        }

        public double CalculateDTWDistance(Timeseries seriesTimeWindow1, Timeseries seriesTimeWindow2)
        {
            List<Timeseries> euclidean_Distances = new List<Timeseries>();
            List<Timeseries> DTW_Distances = new List<Timeseries>();
            for (int n = 0; n < seriesTimeWindow1.length; n++)
            {
                List<double> yDistances = new List<double>();
                List<double> dtwCol = new List<double>();
                for (int m = 0; m < seriesTimeWindow2.length; m++)
                {
                    yDistances.Add(Math.Pow((seriesTimeWindow1.yValues[n] - seriesTimeWindow2.yValues[m]), 2));
                    dtwCol.Add(0.0);
                }
                euclidean_Distances.Add(new Timeseries(seriesTimeWindow1.dateTimes, yDistances));
                DTW_Distances.Add(new Timeseries(seriesTimeWindow1.dateTimes, dtwCol));
            }

            for (int m = 0; m < seriesTimeWindow2.length; m++)
            {
                List<double> dtwCol = new List<double>();
                for (int n = 0; n < seriesTimeWindow1.length; n++)
                {
                    List<double> mins = new List<double>();
                    if ((n > 0) && (m > 0))
                    {
                        mins.Add(DTW_Distances[n - 1].yValues[m]);
                        mins.Add(DTW_Distances[n - 1].yValues[m - 1]);
                        mins.Add(DTW_Distances[n].yValues[m - 1]);
                        DTW_Distances[n].yValues[m] = euclidean_Distances[n].yValues[m] + mins.Min();
                    }
                    else if ((n > 0) && (m == 0))
                    {
                        DTW_Distances[n].yValues[0] = euclidean_Distances[n].yValues[0] + DTW_Distances[n - 1].yValues[0];
                    }
                    else if ((n == 0) && (m > 0))
                    {
                        DTW_Distances[0].yValues[m] = euclidean_Distances[0].yValues[m] + DTW_Distances[0].yValues[m - 1];
                    }
                    else if ((m == 0) && (n == 0))
                    {
                        DTW_Distances[0].yValues[0] = euclidean_Distances[0].yValues[0];
                    }
                }
            }

            int N = seriesTimeWindow1.length - 1;
            int M = seriesTimeWindow2.length - 1;
            int minIndex = 0;
            List<Pair> warpPath = new List<Pair>();
            while ((N + M) > 0)
            {
                if (N - 1 == 0)
                {
                    M = M - 1;
                }
                else
                {
                    if (M - 1 == 0)
                    {
                        N = N - 1;
                    }
                    else
                    {
                        List<double> mins = new List<double>();
                        mins.Add(DTW_Distances[N - 1].yValues[M]);
                        mins.Add(DTW_Distances[N].yValues[M - 1]);
                        mins.Add(DTW_Distances[N - 1].yValues[M - 1]);
                        minIndex = mins.IndexOf(mins.Min());
                    }

                    switch (minIndex)
                    {
                        case 0:
                            N--;
                            break;
                        case 1:
                            M--;
                            break;
                        case 2:
                            N--;
                            M--;
                            break;
                    }
                }

                if ((N <= 0) || (M <= 0)) break;
                Pair p = new Pair(N, M);
                warpPath.Add(p);
            }

            double DTWdistance = 0;
            double warpPathLength = 0;
            double diagonalLength = (warpPath.Count - 1)* Math.Sqrt(2); // length of the diagonal path
            for (int k = 0; k < warpPath.Count - 1; k++)
            {
                warpPathLength += Math.Sqrt(Math.Pow((warpPath[k].value1 - warpPath[k + 1].value1), 2) + Math.Pow((warpPath[k].value2 - warpPath[k + 1].value2), 2)); // warpPath length
            }
            DTWdistance = Math.Abs(warpPathLength - diagonalLength);            

            return DTWdistance;
        }

        private void CreateClusters(Location location, Parameters parameters)
        {
            // use the alglib library to cluster sample timeseries based on the Distance Matrix
            alglib.clusterizerstate s;
            alglib.ahcreport report;

            alglib.clusterizercreate(out s);
            double[,] distanceMatrix = new double[DTWdistancesEvolution[location.name].Count, DTWdistancesEvolution[location.name].Count];
            for (int i = 0; i < DTWdistancesEvolution[location.name].Count; i++)
            {
                for (int j = 0; j < DTWdistancesEvolution[location.name][i].Count; j++)
                {
                    distanceMatrix[i, j] = DTWdistancesEvolution[location.name][i][j];
                }
            }

            DataManager.WriteDistanceMatrix(location, parameters, DTWdistancesEvolution, timeSamples, timeWindowNumber);
            if (timeWindowNumber >= parameters.trainTimeWindowsNumber)
            {
                StoreEvolution(location, parameters);
            }

            alglib.clusterizersetdistances(s, distanceMatrix, true);
            alglib.clusterizersetahcalgo(s, 1);
            alglib.clusterizerrunahc(s, out report);
        }

        public void StoreEvolution(Location location, Parameters parameters)
        {
            List<double> distancesPreviousTimeWindows = new List<double>();
            List<double> distancesNewTimeWindow = new List<double>();
            int distancesEvolutionSize = DTWdistancesEvolution[location.name].Count;
            for (int i = 0; i < distancesEvolutionSize; i++)
            {
                for (int j = i + 1; j < distancesEvolutionSize; j++)
                {
                    if (j != distancesEvolutionSize - 1)
                    {
                        distancesPreviousTimeWindows.Add(DTWdistancesEvolution[location.name][i][j]);
                    }
                    else
                    {
                        distancesNewTimeWindow.Add(DTWdistancesEvolution[location.name][i][j]);
                    }
                }
            }

            // calculate the median of the distances up to the previous timewindow
            double medianPreviousTimewindows = CalculateMedian(distancesPreviousTimeWindows.ToList());

            // calculate the absolute deviations associated with the latest timewindow
            List<double> absoluteDeviations = new List<double>();
            for (int i = 0; i < distancesNewTimeWindow.Count; i++)
            {
                absoluteDeviations.Add(Math.Abs(distancesNewTimeWindow[i] - medianPreviousTimewindows));
            }
            minAbsoluteDeviations.Add(absoluteDeviations.Min());

            DetectOutliers(minAbsoluteDeviations, parameters);

            List<double> minAbsoluteDeviationsWithoutOutliers = new List<double>();
            if(outliers.Count > 0)
            {
                for (int i = 0; i < minAbsoluteDeviations.Count; i++)
                {
                    for (int j = 0; j < outliers.Count; j++)
                    {
                        if (i != outliers[j])
                        {
                            minAbsoluteDeviationsWithoutOutliers.Add(minAbsoluteDeviations[i]);
                            break;
                        }
                    }
                }
                double noOutliersMedian = CalculateMedian(minAbsoluteDeviationsWithoutOutliers);

                minAbsoluteDeviationsWithoutOutliers = new List<double>(minAbsoluteDeviations);
                for (int i = 0; i < minAbsoluteDeviations.Count; i++)
                {
                    for (int j = 0; j < outliers.Count; j++)
                    {
                        if (i == outliers[j])
                        {
                            minAbsoluteDeviationsWithoutOutliers[i] = noOutliersMedian;
                            break;
                        }
                    }
                }
            }

            DataManager.WriteTimeEvolution(location, timeWindowNumber, parameters, minAbsoluteDeviations, minAbsoluteDeviationsWithoutOutliers, outliers);
        }

        private Timeseries ShiftTimeseries(Timeseries timeWindow1Timeseries, Timeseries timeWindow2Timeseries)
        {
            int daysDifference = Convert.ToInt32(timeWindow2Timeseries.dateTimes.First().Subtract(timeWindow1Timeseries.dateTimes.First()).TotalDays);

            List<DateTime> timeWindow2TimeseriesDatetimesShifted = new List<DateTime>();
            List<double> timeWindow2TimeseriesValues = new List<double>();
            foreach (DateTime dt in timeWindow2Timeseries.dateTimes)
            {
                timeWindow2TimeseriesDatetimesShifted.Add(dt.AddDays(-daysDifference));
            }
            timeWindow2TimeseriesValues.AddRange(timeWindow2Timeseries.yValues);

            return new Timeseries(timeWindow2TimeseriesDatetimesShifted, timeWindow2TimeseriesValues);
        }

        private List<DateTime> CreateDatetimesAggregate(Timeseries timeseries1, Timeseries timeseries2)
        {
            List<DateTime> dateTimes = new List<DateTime>();
            dateTimes.AddRange(timeseries1.dateTimes);
            dateTimes.AddRange(timeseries2.dateTimes);

            List<DateTime> distinctDateTimes = dateTimes.Distinct().ToList();
            distinctDateTimes.Sort((a, b) => a.CompareTo(b));

            return distinctDateTimes;
        }

        private Timeseries StretchTimeseries(Timeseries timeWindowTimeseries, List<DateTime> aggregateDatetimes)
        {
            List<DateTime> dateTimes = new List<DateTime>();
            List<double> values = new List<double>();
            int idx = 0;

            foreach (DateTime dt in aggregateDatetimes)
            {
                dateTimes.Add(dt);
                if (timeWindowTimeseries.dateTimes.Contains(dt))
                {
                    values.Add(timeWindowTimeseries.yValues[idx]);
                    idx++;
                }
                else
                {
                    double value1 = 0;
                    if (idx - 1 < 0)
                    {
                        value1 = timeWindowTimeseries.yValues[idx];
                    }
                    else
                    {
                        value1 = timeWindowTimeseries.yValues[idx - 1];
                    }
                    double value2 = 0;
                    if (idx >= timeWindowTimeseries.yValues.Count)
                    {
                        value2 = timeWindowTimeseries.yValues[idx - 1];                        
                    }
                    else
                    {
                        value2 = timeWindowTimeseries.yValues[idx];
                    }

                    double value = (value1 + value2) / 2;
                    values.Add(value);
                }
            }

            return new Timeseries(dateTimes, values);
        }

        private void DetectOutliers(List<double> values, Parameters parameters)
        {
            double median = CalculateMedian(values);

            List<double> absoluteDeviations = new List<double>();
            for (int i = 0; i < values.Count; i++)
            {
                absoluteDeviations.Add(Math.Abs(values[i] - median));
            }

            double MAD = CalculateMedian(absoluteDeviations);

            outliers = new List<int>();
            if(MAD > 0)
            {
                List<double> distancesFromMedian = new List<double>();
                for (int i = 0; i < absoluteDeviations.Count; i++)
                {
                    if (values[i] == median)
                    {
                        distancesFromMedian.Add(0.0);
                    }
                    else
                    {
                        distancesFromMedian.Add(absoluteDeviations[i] / MAD);
                    }
                }

                for (int i = 0; i < distancesFromMedian.Count; i++)
                {
                    if (distancesFromMedian[i] > parameters.anomalyThreshold)
                    {
                        outliers.Add(i);
                    }
                }
            }
        }

        private double CalculateMedian(List<double> values)
        {
            List<double> sortedValues = values.ToList();
            sortedValues.Sort();

            double median;

            int halfIndex = sortedValues.Count / 2;
            if ((sortedValues.Count % 2) == 0)
            {
                median = ((sortedValues.ElementAt(halfIndex) + sortedValues.ElementAt((halfIndex - 1))) / 2);
            }
            else
            {
                median = sortedValues.ElementAt(halfIndex);
            }

            return median;
        }

        public int timeWindowNumber;
        private List<double> minAbsoluteDeviations = new List<double>();
        private List<int> outliers = new List<int>();
        private List<Location> locations = new List<Location>();
        private Dictionary<string, List<List<double>>> DTWdistancesEvolution = new Dictionary<string, List<List<double>>>();
        private Dictionary<string, List<Timeseries>> timeSamples = new Dictionary<string, List<Timeseries>>();
    }
}

