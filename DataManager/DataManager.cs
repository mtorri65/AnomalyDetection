using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.IO.Compression;
using DataStructures;
using ConfigurationParameters;

namespace DataManagement
{
    public class DataManager
    {
        public static Timeseries ReadTimeSeries(string locationName, int timeWindowNumber)
        {
            string timeseriesFilePath;
            timeseriesFilePath = @".\Output\Support\TimeWindow" + timeWindowNumber + @"\" + locationName + @"\" + timeWindowNumber + ".csv";

            Decompress(new FileInfo(timeseriesFilePath));

            string input;
            List<DateTime> dateTimes = new List<DateTime>();
            List<double> values = new List<double>();
            using (StreamReader sr = new StreamReader(timeseriesFilePath))
            {
                while ((input = sr.ReadLine()) != null)
                {
                    string line = input.TrimEnd(',');
                    string[] tokens = line.Split(',');

                    string[] tokensDateTime = tokens[0].Split(' ');
                    string[] tokensDate = tokensDateTime[0].Split('/');
                    int YYYY = Convert.ToInt32(tokensDate[0]);
                    int MM = Convert.ToInt32(tokensDate[1]);
                    int DD = Convert.ToInt32(tokensDate[2]);

                    string[] tokensTime = tokensDateTime[1].Split(':');
                    int hh = Convert.ToInt32(tokensTime[0]);
                    int mm = Convert.ToInt32(tokensTime[1]);
                    int ss = Convert.ToInt32(tokensTime[2]);

                    double value = Convert.ToDouble(tokens[1]);

                    dateTimes.Add(new DateTime(YYYY, MM, DD, hh, mm, ss));
                    values.Add(value);
                }
            }

            Compress(new FileInfo(timeseriesFilePath));

            return new Timeseries(dateTimes, values);
        }

        public static void CleanUpOutput(string folder)
        {
            if (Directory.Exists(folder))
            {
                DirectoryInfo dirInfo = new DirectoryInfo(folder);
                DirectoryInfo[] folderList = dirInfo.GetDirectories();
                foreach (var subDir in folderList)
                {
                    try
                    {
                        subDir.Delete(true);
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine(e.StackTrace);
                        continue;
                    }
                }

                IEnumerable<FileInfo> files = dirInfo.EnumerateFiles();
                if (files.Count() > 0)
                {
                    foreach (FileInfo f in files)
                    {
                        try
                        {
                            f.Delete();
                        }
                        catch (Exception e)
                        {
                            Console.WriteLine(e.StackTrace);
                            continue;
                        }
                    }
                }
            }
        }

        public static void CompressZediAccessData()
        {
            string[] fileEntries = Directory.GetFiles(@".\ZediAccessData");
            FileInfo fileToCompress = null;
            foreach (string fileName in fileEntries)
            {
                fileToCompress = new FileInfo(fileName);
                if (fileToCompress.Extension == ".csv")
                {
                    Compress(fileToCompress);
                }
            }
        }

        public static void CreateHistorianDatabase(List<string> locations)
        {
            if(!Directory.Exists(@".\HistorianDatabase\"))
            {
                Directory.CreateDirectory(@".\HistorianDatabase\");
            }

            string[] filePaths = Directory.GetFiles(@".\HistorianDatabase\");
            foreach (string filePath in filePaths)
            {
                File.Delete(filePath);
            }

            foreach (string location in locations)
            {
                using (FileStream outFile = File.Create(@".\HistorianDatabase\" + location))
                {
                }
                Compress(new FileInfo(@".\HistorianDatabase\" + location));
            }
        }

        public static bool RetrieveLastTimewindowData(string locationName, int timeWindowNumber, int physicalQuantityToMonitor, Parameters parameters)
        {
            bool isThereData = true;

            List<KeyValuePair<DateTime, List<double>>> data = InjectSensorDataInHistorianDatabase(locationName, timeWindowNumber);
            //            if (data.Count == 0)
            if (data.Count < .90 * parameters.timeWindowLength / parameters.reportingFrequency)
            {
//                string fileName = @".\Output\Support\Train\" + locationName + "_log.csv";
//                using (StreamWriter sw = new StreamWriter(fileName, true))
//                {
//                    sw.WriteLine("NO DATA");
                    isThereData = false;
//                }
            }
            else
            {
                if (File.Exists(@".\HistorianDatabase\" + locationName + ".gz"))
                {
                    Decompress(new FileInfo(@".\HistorianDatabase\" + locationName));
                }

                using (StreamWriter sw = new StreamWriter(@".\HistorianDatabase\" + locationName))
                {
                    foreach (KeyValuePair<DateTime, List<double>> kwp in data)
                    {
                        string dateTime = kwp.Key.ToString("yyyy/MM/dd HH:mm:ss").Replace("-", "/") + ",";
                        string values = "";
                        foreach (double d in kwp.Value)
                        {
                            values += d.ToString() + ",";
                        }
                        string record = dateTime + values;

                        sw.WriteLine(record);
                    }
                }

                if (File.Exists(@".\HistorianDatabase\" + locationName))
                {
                    Compress(new FileInfo(@".\HistorianDatabase\" + locationName));
                }
                CreateLocationSupportFiles(data, locationName, timeWindowNumber, physicalQuantityToMonitor, parameters);
            }

            return isThereData;
        }

        public static void WriteDistanceMatrix(Location location, Parameters parameters, 
                                                Dictionary<string, List<List<double>>> distancesEvolution,
                                                Dictionary<string, List<Timeseries>> trainSamples, int timeWindowNumber)
        {
            using (StreamWriter sw = new StreamWriter(@".\Output\Analysis\" + location.name + @"\DistanceTimeMatrix" + timeWindowNumber + ".csv", false))
            {
                for (int row = -1; row < trainSamples[location.name].Count; row++)
                {   // row = -1 is the header
                    string entry = "";
                    for (int col = 0; col < trainSamples[location.name].Count; col++)
                    {
                        string end = ",";
                        if (col == trainSamples[location.name].Count - 1)
                        {
                            end = "\n";
                        }
                        if (row == -1)
                        {
                            entry += col + end;
                        }
                        else
                        {
                            entry += distancesEvolution[location.name][row][col] + end;
                        }
                    }
                    sw.Write(entry);
                }
            }
        }

        public static void WriteTimeEvolution(Location location, int timeWindowNumber, Parameters parameters, 
                                              List<double> minAbsoluteDeviations, List<double> minAbsoluteDeviationsWithoutOutliers, List<int> outliers)
        {
            using (StreamWriter sw = new StreamWriter(@".\Output\Analysis\" + location.name + @"\TimeEvolution_" + location.name, true))
            {
                string outlier = minAbsoluteDeviations.Last().ToString();
                if (outliers.Count > 0)
                {
                    outlier = minAbsoluteDeviationsWithoutOutliers.Last() + ",";
                }
                sw.WriteLine(timeWindowNumber + "," + minAbsoluteDeviations.Last() + "," + outlier);
            }
        }

        private static void CreateLocationSupportFiles(List<KeyValuePair<DateTime, List<double>>> data,
                                                        string location,                                                
                                                        int timeWindowNumber,
                                                        int physicalQuantityToMonitor,
                                                        Parameters parameters)
        {
            string locationSupportFolder = @".\Output\Support\Timewindow" + timeWindowNumber + @"\" + location;
            Directory.CreateDirectory(locationSupportFolder);
            string locationAnalysisFolder = @".\Output\Analysis\" + location;
            Directory.CreateDirectory(locationAnalysisFolder);

            DateTime firstTimestamp = data.First().Key;
            DateTime lastTimestamp = data.Last().Key;
            List<DateTime> endOfTimeWindows = new List<DateTime>();
            List<StreamWriter> swTimeseriesCSVs = new List<StreamWriter>();

            int timeWindowLengthInDays = parameters.timeWindowLength / 86400;
            TimeSpan numberOfDaysFromStart = new TimeSpan(timeWindowLengthInDays, 0, 0, 0);
            endOfTimeWindows.Add(lastTimestamp);
            swTimeseriesCSVs.Add(new StreamWriter(locationSupportFolder + @"\" + timeWindowNumber + ".csv"));

            for (int i = 0; i < data.Count; i++)
            {
                int timeWindowLength = 0;
                DateTime dt = data[i].Key;
                for(int w = 0; w < endOfTimeWindows.Count - 1; w++)
                {
                    if(dt <= endOfTimeWindows[0])
                    {
                        timeWindowLength = 0;
                        break;
                    }
                    else if((endOfTimeWindows[w] <= dt) && (dt <= endOfTimeWindows[w + 1]))
                    {
                        timeWindowLength = w + 1;
                        break;
                    }
                }

                string sdt = dt.ToString("yyyy/MM/dd HH:mm:ss").Replace("-", "/");
                double dp = data[i].Value[0];
                double ft = data[i].Value[1];
                double gf = data[i].Value[2];
                double sp = data[i].Value[3];

                double analyzedQuantity = data[i].Value[physicalQuantityToMonitor - 1];

                swTimeseriesCSVs[timeWindowLength].WriteLine(sdt + "," + analyzedQuantity + ",");
            }

            foreach(StreamWriter sw in swTimeseriesCSVs)
            {
                sw.Close();
            }

            for(int m = 1; m <= endOfTimeWindows.Count; m++)
            {
                CompressFiles(locationSupportFolder, @"\" + timeWindowNumber, ".csv");
            }
        }

        private static void Compress(FileInfo fileToCompress)
        {
            using (FileStream originalFileStream = fileToCompress.OpenRead())
            {
                string currentFileName = fileToCompress.FullName;
                string newFileName = fileToCompress.FullName + ".gz";

                using (FileStream compressedFileStream = File.Create(newFileName))
                {
                    using (GZipStream compressionStream = new GZipStream(compressedFileStream, CompressionMode.Compress))
                    {
                        originalFileStream.CopyTo(compressionStream);
                    }
                }
            }
            try
            {
                File.Delete(fileToCompress.FullName);
            }
            catch (Exception e)
            {
                Console.WriteLine(e.StackTrace);
                return;
            }
        }

        private static void Decompress(FileInfo fileToDecompress)
        {
            FileInfo fileToDecompressWithExtension = new FileInfo(fileToDecompress.FullName + ".gz");
            using (FileStream originalFileStream = fileToDecompressWithExtension.OpenRead())
            {
                string currentFileName = fileToDecompressWithExtension.FullName;
                string newFileName = currentFileName.Remove(currentFileName.Length - fileToDecompressWithExtension.Extension.Length);

                using (FileStream decompressedFileStream = File.Create(newFileName))
                {
                    using (GZipStream decompressionStream = new GZipStream(originalFileStream, CompressionMode.Decompress))
                    {
                        decompressionStream.CopyTo(decompressedFileStream);
                    }
                }
            }
            try
            {
                File.Delete(fileToDecompressWithExtension.FullName);
            }
            catch (Exception e)
            {
                Console.WriteLine(e.StackTrace);
                return;
            }
        }

        private static void CompressFiles(string folder, string fileName, string extension)
        {
            FileInfo fInfo = new FileInfo(folder + fileName + extension);
            Compress(fInfo);
        }

        private static List<KeyValuePair<DateTime, List<double>>> InjectSensorDataInHistorianDatabase(string location, int timeWindowNumber)
        {
            string timeSeries = @".\ZediAccessData\" + location;

            if (File.Exists(timeSeries + ".gz"))
            {
                Decompress(new FileInfo(timeSeries));
            }

            int timeseriesLength = File.ReadAllLines(timeSeries).Length;
            if (timeseriesLength <= 2)
            {
                Console.WriteLine("File " + timeSeries + " does not contain any timestamped value!");
            }

            List<KeyValuePair<DateTime, List<double>>> data = null;
            Console.WriteLine("Day # = " + timeWindowNumber);
            Console.WriteLine("Location = " + location);
            DateTime firstTimestamp = new DateTime();
            using (StreamReader sr = new StreamReader(timeSeries))
            {
                string input;
                int n = 0;
                int offset = 2;
                data = new List<KeyValuePair<DateTime, List<double>>>();

                while ((input = sr.ReadLine()) != null)
                {
                    ++n;
                    if (n > offset)
                    {
                        string[] tokens = input.Split(',');
                        /*                                int YYYY = Convert.ToInt32(tokens[0].Substring(0, 4));
                                                    int MM = Convert.ToInt32(tokens[0].Substring(4, 2));
                                                    int DD = Convert.ToInt32(tokens[0].Substring(6, 2));
                                                    string time = tokens[0].Substring(9); */
                        string[] tokensDateTime = tokens[0].Split(' ');
                        string[] tokensDate = tokensDateTime[0].Split('/');
                        int MM = Convert.ToInt32(tokensDate[0]);
                        int DD = Convert.ToInt32(tokensDate[1]);
                        int YYYY = Convert.ToInt32(tokensDate[2]);

                        string[] tokensTime = tokensDateTime[1].Split(':');
                        int hh = Convert.ToInt32(tokensTime[0]);
                        int mm = Convert.ToInt32(tokensTime[1]);
                        int ss = 0;
                        if (tokensTime.Count() == 3)
                        {
                            ss = Convert.ToInt32(tokensTime[2]);
                        }

                        DateTime timestamp = new DateTime(YYYY, MM, DD, hh, mm, ss);
                        TimeSpan firstTimeWindowTimestamp = new TimeSpan(0, 0, 0, 0);
                        if ((n - 1) == offset)
                        {
                            firstTimestamp = timestamp + firstTimeWindowTimestamp;
                        }
                        TimeSpan timePassedFromTimeWindowFirstTimeStamp = timestamp.Subtract(firstTimestamp);
                        TimeSpan daysBeforeThisTimeWindowInLocationData = new TimeSpan(timeWindowNumber - 1, 0, 0, 0);
                        if (TimeSpan.Compare(timePassedFromTimeWindowFirstTimeStamp, daysBeforeThisTimeWindowInLocationData) < 0)
                        {
                            continue;
                        }
                        TimeSpan daysAfterThisTimeWindowInLocationData = new TimeSpan(timeWindowNumber, 0, 0, 0);
                        if (TimeSpan.Compare(timePassedFromTimeWindowFirstTimeStamp, daysAfterThisTimeWindowInLocationData) >= 0)
                        {
                            break;
                        }

                        List<double> values = new List<double>();
                        for (int i = 1; i < tokens.Length; i++)
                        {
                            values.Add(Convert.ToDouble(tokens[i]));
                        }

                        data.Add(new KeyValuePair<DateTime, List<double>>(timestamp, values));
                    }
                }
            }

            if (File.Exists(timeSeries))
            {
                Compress(new FileInfo(timeSeries));
            }

            return data;
        }
    }
}
