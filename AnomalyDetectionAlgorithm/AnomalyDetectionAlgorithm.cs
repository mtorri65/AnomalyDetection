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
            this.frequencyWindowNumber = timeWindowNumber;
            
            foreach(string locationName in locationNames)
            {
                Location location = new Location(locationName);
                locations.Add(location);
            }
        }

        public void Train_Time(Parameters parameters, string locationName)
        {
            Location location = locations.Find(x => x.name == locationName);
            timeSamples[location.name] = new List<Timeseries>();
            for (int timeWindowNumber = 1; timeWindowNumber <= parameters.trainTimeWindowsNumber; timeWindowNumber++)
            {
                Series series = DataManager.ReadTimeSeries(location.name, timeWindowNumber/*, "TRAIN"*/);
                timeSamples[location.name].Add(series.toTimeseries());
            }

            // calculate the DTW distance between each pair of timeseries in the training set (Distance Matrix)
            List<List<double>> DTWtimeWindowDistances = new List<List<double>>();
            for (int i = 0; i < timeSamples[location.name].Count; i++)
            {
                DTWtimeWindowDistances.Add(new List<double>(new double[timeSamples[location.name].Count]));
            }
            DTWdistancesEvolution_Time.Add(location.name, DTWtimeWindowDistances);

            if (parameters.distMatrixExists == "yes")
            {
//                DataManager.ReadDistanceMatrix(location, DTWdistancesEvolution_Time, timeWindowNumber);
            }
            else
            {
                for (int i = 0; i < timeSamples[location.name].Count; i++)
                {
                    for (int j = i; j < timeSamples[location.name].Count; j++)
                    {
                        if (j != i  )
                        {
                            Timeseries sampleShifted = ShiftTimeseries(timeSamples[location.name][i], timeSamples[location.name][j]);
                            List<DateTime> aggregateDatetimes = CreateDatetimesAggregate(timeSamples[location.name][i], sampleShifted);
                            Timeseries sampleStretched1 = StretchTimeseries(timeSamples[location.name][i], aggregateDatetimes);
                            Timeseries sampleStretched2 = StretchTimeseries(sampleShifted, aggregateDatetimes);

                            DTWdistancesEvolution_Time[location.name][i][j] = CalculateDTWDistance_Time(sampleStretched1, sampleStretched2);
                            DTWdistancesEvolution_Time[location.name][j][i] = DTWdistancesEvolution_Time[location.name][i][j];
                        }
                        else
                        {
                            DTWdistancesEvolution_Time[location.name][i][j] = 0.0;
                        }
                    }
                } 
            }

            CreateClusters_Time(location, parameters);
        }

/*        public void Train_Frequency(Parameters parameters)
        {
            foreach (Location location in locations)
            {
                frequencySamples[location.name] = new List<FrequencySeries>();
                for (int timeWindowNumber = 1; timeWindowNumber <= parameters.trainTimeWindowsNumber; timeWindowNumber++)
                {
                    Series series = DataManager.ReadTimeSeries(location.name, timeWindowNumber);
                    frequencySamples[location.name].Add(series.toFrequencySeries());
                }

                // calculate the DTW distance between each pair of timeseries in the training set (Distance Matrix)
                List<List<double>> DTWfrequencyWindowDistances = new List<List<double>>();
                for (int i = 0; i < frequencySamples[location.name].Count; i++)
                {
                    DTWfrequencyWindowDistances.Add(new List<double>(new double[frequencySamples[location.name].Count]));
                }
                List<List<List<double>>> DTWfrequencyDistances = new List<List<List<double>>>();

                DTWfrequencyDistances.Add(DTWfrequencyWindowDistances);
                DTWdistancesEvolution_Frequency.Add(location.name, DTWfrequencyDistances);

                if (parameters.distMatrixExists == "yes")
                {
                    DataManager.ReadDistanceMatrix(location, DTWdistancesEvolution_Time, timeWindowNumber);
                }
                else
                {
                    for (int i = 0; i < timeSamples[location.name].Count; i++)
                    {
                        for (int j = i; j < timeSamples[location.name].Count; j++)
                        {
                            if (j != i)
                            {
                                DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][i][j] = CalculateDTWDistance_Frequency(frequencySamples[location.name][i], frequencySamples[location.name][j]);
                                DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][j][i] = DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][i][j];
                            }
                            else
                            {
                                DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][i][j] = 0.0;
                                DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][i][j] = 0.0;
                            }
                        }
                    }
                    for (int i = 0; i < timeSamples[location.name].Count; i++)
                    {
                        for (int j = i; j < timeSamples[location.name].Count; j++)
                        {
                            if (j != i)
                            {
                                DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][i][j] = CalculateDTWDistance_Frequency(frequencySamples[location.name][i], frequencySamples[location.name][j]);
                                DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][j][i] = DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][i][j];
                            }
                            else
                            {
                                DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][i][j] = 0.0;
                                DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - 1][i][j] = 0.0;
                            }
                        }
                    }
                }

                CreateClusters_Frequency(location, parameters);
            }
        }*/

        private void CreateClusters_Time(Location location, Parameters parameters)
        {
            // use the alglib library to cluster sample timeseries based on the Distance Matrix
            alglib.clusterizerstate s;
            alglib.ahcreport report;

            alglib.clusterizercreate(out s);
            double[,] distanceMatrix = new double[DTWdistancesEvolution_Time[location.name].Count, DTWdistancesEvolution_Time[location.name].Count];
            for (int i = 0; i < DTWdistancesEvolution_Time[location.name].Count; i++)
            {
                for (int j = 0; j < DTWdistancesEvolution_Time[location.name][i].Count; j++)
                {
                    distanceMatrix[i, j] = DTWdistancesEvolution_Time[location.name][i][j];
                }
            }

            DataManager.WriteDistanceMatrix_Time(location, parameters, DTWdistancesEvolution_Time, timeSamples, timeWindowNumber);
            if(timeWindowNumber >= parameters.trainTimeWindowsNumber)
            {
                StoreEvolution_Time(location, parameters);
            }

            alglib.clusterizersetdistances(s, distanceMatrix, true);
            alglib.clusterizersetahcalgo(s, 1);
            alglib.clusterizerrunahc(s, out report);
        }

/*        private void CreateClusters_Frequency(Location location, Parameters parameters)
        {
            // use the alglib library to cluster sample timeseries based on the Distance Matrix
            alglib.clusterizerstate s;
            alglib.ahcreport report;

            alglib.clusterizercreate(out s);
            int offset = 0;
            if (parameters.train == "yes")
            {
                offset = 1;
            }
            double[,] distanceMatrix = new double[DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber - offset].Count, DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber - offset].Count];
            for (int i = 0; i < DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber - offset].Count; i++)
            {
                for (int j = 0; j < DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber - offset][i].Count; j++)
                {
                    distanceMatrix[i, j] = DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber - offset][i][j];
                }
            }

            DataManager.WriteDistanceMatrix_Frequency(location, parameters, offset, DTWdistancesEvolution_Frequency, frequencySamples, frequencyWindowNumber);

            StoreEvolution_Frequency(location, frequencyWindowNumber, parameters, offset);

            alglib.clusterizersetdistances(s, distanceMatrix, true);
            alglib.clusterizersetahcalgo(s, 1);
            alglib.clusterizerrunahc(s, out report);
        }*/

        public void Detect_Time(Parameters parameters, string locationName)
        {
            Location location = locations.Find(x => x.name == locationName);

            Timeseries thisTimeWindowSample = DataManager.ReadTimeSeries(location.name, timeWindowNumber/*, "TEST"*/).toTimeseries();
            if(timeWindowNumber == 1)
            {
                List<List<double>> DTWtimeWindowDistances = new List<List<double>>();
                DTWtimeWindowDistances.Add(new List<double>(new double[1]));
                DTWdistancesEvolution_Time.Add(location.name, DTWtimeWindowDistances);

                List<Timeseries> lts = new List<Timeseries>();
                lts.Add(thisTimeWindowSample);
                timeSamples.Add(location.name, lts);

                return;
            }

            // add one extra column (and one extra row) to the Distance Matrix by including this timeWindowNumber's timeseries in it
            for (int row = 0; row < DTWdistancesEvolution_Time[location.name].Count; row++)
            {
                DTWdistancesEvolution_Time[location.name][row].Add(0.0);
            }
            DTWdistancesEvolution_Time[location.name].Add((new double[DTWdistancesEvolution_Time[location.name].Count + 1]).ToList());

            for (int i = 0; i < DTWdistancesEvolution_Time[location.name].Count; i++)
            {
                if (i < DTWdistancesEvolution_Time[location.name].Count - 1)
                {
                    Timeseries sampleShifted = ShiftTimeseries(timeSamples[location.name][i], thisTimeWindowSample);
                    List<DateTime> aggregateDatetimes = CreateDatetimesAggregate(timeSamples[location.name][i], sampleShifted);
                    Timeseries sampleStretched1 = StretchTimeseries(timeSamples[location.name][i], aggregateDatetimes);
                    Timeseries sampleStretched2 = StretchTimeseries(sampleShifted, aggregateDatetimes);

                    DTWdistancesEvolution_Time[location.name][i][DTWdistancesEvolution_Time[location.name].Count - 1] = CalculateDTWDistance_Time(sampleStretched1, sampleStretched2);
                    DTWdistancesEvolution_Time[location.name][DTWdistancesEvolution_Time[location.name].Count - 1][i] = DTWdistancesEvolution_Time[location.name][i][DTWdistancesEvolution_Time[location.name].Count - 1];
                }
                else
                {
                    DTWdistancesEvolution_Time[location.name][i][DTWdistancesEvolution_Time[location.name].Count - 1] = 0.0; // diagonal term
                }
            }
            timeSamples[location.name].Add(thisTimeWindowSample);
            CreateClusters_Time(location, parameters);
        }

/*        public bool Detect_Frequency(Parameters parameters)
        {
            bool successfulCompletion = false;

            foreach (Location location in locations)
            {
                // for the new frequencyWindowNumber, add the distanceMatrix from last frequencyWindowNumber
                DTWdistancesEvolution_Frequency[location.name].Add(new List<List<double>>(DTWdistancesEvolution_Frequency[location.name][DTWdistancesEvolution_Frequency[location.name].Count - 1]));
                FrequencySeries thisFrequencyWindowSample = DataManager.ReadTimeSeries(location.name, frequencyWindowNumber).toFrequencySeries();

                // add one extra column (and one extra row) to the Distance Matrix by including this frequencyWindowNumber's timeseries in it
                for (int row = 0; row < /DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Count; row++)
                {
                    DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber][row].Add(0.0);
                }
                DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Add((new double[DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Count + 1]).ToList());

                for (int i = 0; i < DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Count; i++)
                {
                    if (i < DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Count - 1)
                    {
                        DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber][i][DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Count - 1] = CalculateDTWDistance_Frequency(frequencySamples[location.name][i], thisFrequencyWindowSample);
                        DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber][DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Count - 1][i] = DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber][i][DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Count - 1];
                    }
                    else
                    {
                        DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber][i][DTWdistancesEvolution_Frequency[location.name][frequencyWindowNumber].Count - 1] = 0.0; // diagonal term
                    }
                }
                frequencySamples[location.name].Add(thisFrequencyWindowSample);
                CreateClusters_Frequency(location, parameters);
            }

            ++frequencyWindowNumber;

            successfulCompletion = true;
            return successfulCompletion;
        }*/

        public double CalculateDTWDistance_Time(Timeseries seriesTimeWindow1, Timeseries seriesTimeWindow2)
        {
            List<Timeseries> euclidean_Distances = new List<Timeseries>();
            List<Timeseries> DTW_Distances = new List<Timeseries>();
            for (int n = 0; n < seriesTimeWindow1.length; n++)
            {
                List<double> yDistances = new List<double>();
                List<double> dtwCol = new List<double>();
                for (int m = 0; m < seriesTimeWindow2.length; m++)
                {
                    yDistances.Add(Math.Pow((seriesTimeWindow1.normalized_yValues[n] - seriesTimeWindow2.normalized_yValues[m]), 2));
                    dtwCol.Add(0.0);
                }
                Series s = new Series(seriesTimeWindow1.dateTimes, yDistances);
                Timeseries t = s.toTimeseries();
                euclidean_Distances.Add((new Series(seriesTimeWindow1.dateTimes, yDistances)).toTimeseries());
                DTW_Distances.Add((new Series(seriesTimeWindow1.dateTimes, dtwCol)).toTimeseries());
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

        public double CalculateDTWDistance_Frequency(FrequencySeries seriesTimeWindow1, FrequencySeries seriesTimeWindow2)
        {
            List<FrequencySeries> euclidean_Distances = new List<FrequencySeries>();
            List<FrequencySeries> DTW_Distances = new List<FrequencySeries>();
            for (int n = 0; n < seriesTimeWindow1.length; n++)
            {
                List<double> frequencyDistances = new List<double>();
                List<double> dtwCol = new List<double>();
                for (int m = 0; m < seriesTimeWindow2.length; m++)
                {
                    frequencyDistances.Add(Math.Pow((seriesTimeWindow1.fValues[n] - seriesTimeWindow2.fValues[m]), 2));
                    dtwCol.Add(0.0);
                }
                euclidean_Distances.Add(new FrequencySeries(seriesTimeWindow1.xValues, frequencyDistances));
                DTW_Distances.Add(new FrequencySeries(seriesTimeWindow1.xValues, dtwCol));
            }

            for (int m = 0; m < seriesTimeWindow2.length; m++)
            {
                List<double> dtwCol = new List<double>();
                for (int n = 0; n < seriesTimeWindow1.length; n++)
                {
                    List<double> mins = new List<double>();
                    if ((n > 0) && (m > 0))
                    {
                        mins.Add(DTW_Distances[n - 1].fValues[m]);
                        mins.Add(DTW_Distances[n - 1].fValues[m - 1]);
                        mins.Add(DTW_Distances[n].fValues[m - 1]);
                        DTW_Distances[n].fValues[m] = euclidean_Distances[n].fValues[m] + mins.Min();
                    }
                    else if ((n > 0) && (m == 0))
                    {
                        DTW_Distances[n].fValues[0] = euclidean_Distances[n].fValues[0] + DTW_Distances[n - 1].fValues[0];
                    }
                    else if ((n == 0) && (m > 0))
                    {
                        DTW_Distances[0].fValues[m] = euclidean_Distances[0].fValues[m] + DTW_Distances[0].fValues[m - 1];
                    }
                    else if ((m == 0) && (n == 0))
                    {
                        DTW_Distances[0].fValues[0] = euclidean_Distances[0].fValues[0];
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
                        mins.Add(DTW_Distances[N - 1].fValues[M]);
                        mins.Add(DTW_Distances[N].fValues[M - 1]);
                        mins.Add(DTW_Distances[N - 1].fValues[M - 1]);
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
            double diagonalLength = (warpPath.Count - 1) * Math.Sqrt(2); // length of the diagonal path
            for (int k = 0; k < warpPath.Count - 1; k++)
            {
                warpPathLength += Math.Sqrt(Math.Pow((warpPath[k].value1 - warpPath[k + 1].value1), 2) + Math.Pow((warpPath[k].value2 - warpPath[k + 1].value2), 2)); // warpPath length
            }
            DTWdistance = Math.Abs(warpPathLength - diagonalLength);

            return DTWdistance;
        }

/*        public void DetectAnomaly(Location location, alglib.ahcreport report, Parameters parameters)
        {
            int optimalNumCluster = 2;

            // create the dendrogram
            int[] clustersIdxs;
            int[] cz;
            alglib.clusterizergetkclusters(report, optimalNumCluster, out clustersIdxs, out cz);

            // alglib permutes the timeseries in the training set to build a non-intersecting dendrogram
            int offset = 0; ;
            if (parameters.train == "yes")
            {
                offset = 1;
            }
            // the array trainSamplesPermuted is just a more user-friendly way of referring to the permuted samples (stored in the alglib report.p array)            
            Timeseries[] trainSamplesPermuted = new Timeseries[DTWdistancesEvolution_Time[location.name][timeWindowNumber - offset].Count];
            for (int i = 0; i < timeSamples[location.name].Count; i++)
            {
                trainSamplesPermuted[report.p[i]] = timeSamples[location.name][i];
            }
            
            // fill up the clusters with the samples that belong to each one
            List<List<KeyValuePair<int, Timeseries>>> clusters = new List<List<KeyValuePair<int, Timeseries>>>();
            for (int i = 0; i < optimalNumCluster; i++)
            {
                clusters.Add(new List<KeyValuePair<int, Timeseries>>());
            }
            for (int clusterNumber = 0; clusterNumber < optimalNumCluster; clusterNumber++)
            {
                for (int clusterIdxSampleBelongsTo = 0; clusterIdxSampleBelongsTo < clustersIdxs.Length; clusterIdxSampleBelongsTo++)
                {
                    if (clustersIdxs[clusterIdxSampleBelongsTo] == clusterNumber)
                    {
                        clusters[clusterNumber].Add(new KeyValuePair<int, Timeseries>(clusterIdxSampleBelongsTo, trainSamplesPermuted[clusterIdxSampleBelongsTo]));
                    }
                }
            }

            DataManager.WriteLocationLog(location, timeWindowNumber, report, parameters, timeSamples[location.name], offset, clusters, 
                                            minClusterSize, potentialAnomalyClusterRelativeSize, clusterAnomaliesDistanceRatio, clustersIdxs);
        }*/

        public void StoreEvolution_Time(Location location, /*int timeWindowNumber,*/ Parameters parameters/*, int offset*/)
        {
            List<double> distancesPreviousTimeWindows = new List<double>();
            List<double> distancesNewTimeWindow = new List<double>();
            int distancesEvolutionSize = DTWdistancesEvolution_Time[location.name].Count;
            for (int i = 0; i < distancesEvolutionSize; i++)
            {
                for (int j = i + 1; j < distancesEvolutionSize; j++)
                {
                    if (j != distancesEvolutionSize - 1)
                    {
                        distancesPreviousTimeWindows.Add(DTWdistancesEvolution_Time[location.name][i][j]);
                    }
                    else
                    {
                        distancesNewTimeWindow.Add(DTWdistancesEvolution_Time[location.name][i][j]);
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

            outliers = DetectOutliers(minAbsoluteDeviations, parameters);

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

            DataManager.WriteTimeEvolution_Time(location, timeWindowNumber, parameters, /*offset,*/ minAbsoluteDeviations, minAbsoluteDeviationsWithoutOutliers, outliers);
        }

        public void StoreEvolution_Frequency(Location location, int timeWindowNumber, Parameters parameters, int offset)
        {
            List<double> distancesPreviousTimeWindows = new List<double>();
            List<double> distancesNewTimeWindow = new List<double>();
            int distancesEvolutionSize = DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - offset].Count;
            for (int i = 0; i < distancesEvolutionSize; i++)
            {
                for (int j = i + 1; j < distancesEvolutionSize; j++)
                {
                    if (j != distancesEvolutionSize - 1)
                    {
                        distancesPreviousTimeWindows.Add(DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - offset][i][j]);
                    }
                    else
                    {
                        distancesNewTimeWindow.Add(DTWdistancesEvolution_Frequency[location.name][timeWindowNumber - offset][i][j]);
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

            outliers = DetectOutliers(minAbsoluteDeviations, parameters);

            List<double> minAbsoluteDeviationsWithoutOutliers = new List<double>();
            if (outliers.Count > 0)
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

            DataManager.WriteTimeEvolution_Frequency(location, timeWindowNumber, parameters, offset, minAbsoluteDeviations, minAbsoluteDeviationsWithoutOutliers, outliers);
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

            return new Series(timeWindow2TimeseriesDatetimesShifted, timeWindow2TimeseriesValues).toTimeseries();
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

            return new Series(dateTimes, values).toTimeseries();
        }

        public List<int> DetectOutliers(List<double> values, Parameters parameters)
        {
            double median = CalculateMedian(values);

            List<double> absoluteDeviations = new List<double>();
            for (int i = 0; i < values.Count; i++)
            {
                absoluteDeviations.Add(Math.Abs(values[i] - median));
            }

            double MAD = CalculateMedian(absoluteDeviations);

            List<int> outliers = new List<int>();
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

            return outliers;
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
        public int frequencyWindowNumber;
        public List<double> minAbsoluteDeviations = new List<double>();
        public List<int> outliers = new List<int>();
        private List<Location> locations = new List<Location>();
        private Dictionary<string, List<List<double>>> DTWdistancesEvolution_Time = new Dictionary<string, List<List<double>>>();
        private Dictionary<string, List<List<List<double>>>> DTWdistancesEvolution_Frequency = new Dictionary<string, List<List<List<double>>>>();
        private Dictionary<string, List<Timeseries>> timeSamples = new Dictionary<string, List<Timeseries>>();
        private Dictionary<string, List<FrequencySeries>> frequencySamples = new Dictionary<string, List<FrequencySeries>>();
        private List<int> minClusterSize = new List<int>();
        private List<double> potentialAnomalyClusterRelativeSize = new List<double>();
        private List<double> clusterAnomaliesDistanceRatio = new List<double>();
    }
}

