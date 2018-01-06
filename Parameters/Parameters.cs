using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConfigurationParameters
{
    public class Parameters
    {
        public string restart;
        public string scalingType;
        public int trainTimeWindowsNumber;
        public string distMatrixExists;
        public double anomalyThreshold;
        public double anomaliesClusterSizeThreshold;
        public double mergeDistanceRatioThreshold;
        public string PSD;
    }
}
