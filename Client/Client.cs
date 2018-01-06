using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using DataManagement;
using AnomalyDetection;
using ConfigurationParameters;

namespace Client
{
    class Client
    {
        static void Main(string[] args)
        {
            Parameters parameters = ReadParameters();
            AnomalyDetectionAlgorithm algorithm = null;
            CreateLocations();

            if (parameters.restart.Equals("yes"))
            {
                DataManager.CleanUpOutput(@".\Output\");

                if (locations.Count > 0)
                {
                    DataManager.CreateHistorianDatabase(locations);
                }
                else
                {
                    Console.WriteLine("Locations are not defined!");
                    return;
                }
            }
            algorithm = CreateDetectionAlgorithm();

            foreach (string location in locations)
            {
                bool isThereDataThisTimeWindow = DataManager.RetrieveLastTimewindowData(location, timeWindowNumber, physicalQuantityToMonitor);
                if (isThereDataThisTimeWindow)
                {
                    algorithm.Detect_Time(parameters, location);
                }
            }
            ++algorithm.timeWindowNumber;

            SerializeAlgorithm(algorithm);

            string[] linesInFile = File.ReadAllLines(@".\Input\Parameters.txt");
            linesInFile[0] = "restart\tno";
            File.WriteAllLines(@".\Input\Parameters.txt", linesInFile);
        }

        private static Parameters ReadParameters()
        {
            Parameters parameters = new Parameters();
            string input;
            using (StreamReader sr = new StreamReader(@".\Input\Parameters.txt"))
            {
                while ((input = sr.ReadLine()) != null)
                {
                    int i = input.LastIndexOf('\t');
                    string value = input.Substring(i + 1);
                    if (input.Contains("restart"))
                    {
                        parameters.restart = value;
                        continue;
                    }
                    if (input.Contains("scalingType"))
                    {
                        parameters.scalingType = value;
                        continue;
                    }
                    if (input.Contains("trainTimeWindowsNumber"))
                    {
                        parameters.trainTimeWindowsNumber = Convert.ToInt32(value);
                        continue;
                    }
                    if (input.Contains("distMatrixExists"))
                    {
                        parameters.distMatrixExists = value;
                        continue;
                    }
                    if (input.Contains("anomalyThreshold"))
                    {
                        parameters.anomalyThreshold = Convert.ToDouble(value);
                        continue;
                    }
                    if (input.Contains("anomaliesClusterSizeThreshold"))
                    {
                        parameters.anomaliesClusterSizeThreshold = Convert.ToDouble(value);
                        continue;
                    }
                    if (input.Contains("mergeDistanceRatioThreshold"))
                    {
                        parameters.mergeDistanceRatioThreshold = Convert.ToDouble(value);
                        continue;
                    }
                }
            }
                                
            return parameters;
        }

        private static void CreateLocations()
        {
            string locationNamesFile = @".\Input\Locations.txt";
            using (StreamReader sr = new StreamReader(locationNamesFile))
            {
                string locatioName;
                while ((locatioName = sr.ReadLine()) != null)
                {
                    locations.Add(locatioName);
                }
            }
        }

        private static AnomalyDetectionAlgorithm CreateDetectionAlgorithm()
        {
            AnomalyDetectionAlgorithm algorithm = null;
            if (File.Exists(algorithmSerializedFile))
            {
                algorithm = DeserializeAlgorithm();
                timeWindowNumber = algorithm.timeWindowNumber;
            }
            else
            {
                algorithm = new AnomalyDetectionAlgorithm(locations, timeWindowNumber);
            }

            return algorithm;
        }

        private static void SerializeAlgorithm(AnomalyDetectionAlgorithm algorithm)
        {
            Stream algorithmStream = File.Create(algorithmSerializedFile);
            BinaryFormatter sr = new BinaryFormatter();
            sr.Serialize(algorithmStream, algorithm);
            algorithmStream.Close();
        }

        private static AnomalyDetectionAlgorithm DeserializeAlgorithm()
        {
            Stream algorithmStream = File.OpenRead(algorithmSerializedFile);
            BinaryFormatter deserializer = new BinaryFormatter();
            AnomalyDetectionAlgorithm algorithm = (AnomalyDetectionAlgorithm)deserializer.Deserialize(algorithmStream);
            algorithmStream.Close();
            return algorithm;
        }

        private static List<string> locations = new List<string>();
        private static int physicalQuantityToMonitor = 1;
        private static int timeWindowNumber = 1;

        private const string algorithmSerializedFile = @".\Output\AlgorithmClient.bin";
    }
}
