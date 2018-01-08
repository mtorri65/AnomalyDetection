using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DataStructures
{
    [Serializable]
    public class Timeseries
    {
        public Timeseries()
        {
            dateTimes = new List<DateTime>();
            yValues = new List<double>();
        }

        public Timeseries(List<DateTime> dateTimes, List<double> yValues)
        {
            this.dateTimes = dateTimes;
            this.yValues = yValues;

            this.length = yValues.Count;
        }

        public List<DateTime> dateTimes = new List<DateTime>();
        public List<double> yValues = new List<double>();

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
