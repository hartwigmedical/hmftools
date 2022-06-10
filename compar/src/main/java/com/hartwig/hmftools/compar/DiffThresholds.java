package com.hartwig.hmftools.compar;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.util.Map;

import com.google.common.collect.Maps;

public class DiffThresholds
{
    private final Map<String,ThresholdData> mFieldThresholds;

    // private static final String DEFAULT_FIELD = "DEFAULT";

    private static final String THRESHOLD_ITEM_DELIM = ":";

    public DiffThresholds()
    {
        mFieldThresholds = Maps.newHashMap();
    }

    public boolean hasDifference(final String field, double value1, double value2)
    {
        ThresholdData thresholdData = mFieldThresholds.get(field);

        if(thresholdData == null)
            return false;

        return thresholdData.hasDiff(value1, value2);
    }

    public void addFieldThreshold(final String field, double absoluteDiff, double percentDiff)
    {
        ThresholdType type = absoluteDiff > 0 && percentDiff > 0 ? ThresholdType.ABSOLUTE_AND_PERCENT :
                (absoluteDiff > 0 ? ThresholdType.ABSOLUTE : ThresholdType.PERCENT);

        mFieldThresholds.put(field, new ThresholdData(type, absoluteDiff, percentDiff));
    }

    public void loadConfig(final String configStr)
    {
        if(configStr == null || configStr.isEmpty())
            return;

        // in the form: Field:AbsDiff:PercDiff and separated by ;
        String[] fieldEntries = configStr.split(ITEM_DELIM);

        for(String fieldEntry : fieldEntries)
        {
            String[] thresholdItems = fieldEntry.split(THRESHOLD_ITEM_DELIM, 3);
            String field = thresholdItems[0];
            double absDiff = Double.parseDouble(thresholdItems[1]);
            double percDiff = Double.parseDouble(thresholdItems[2]);
            addFieldThreshold(field, absDiff, percDiff);

            CMP_LOGGER.info("added threshold: field({}) absoluteDiff({}) percentDiff({})", field, absDiff, percDiff);
        }
    }

    private enum ThresholdType
    {
        ABSOLUTE,
        PERCENT,
        ABSOLUTE_AND_PERCENT
    }

    private class ThresholdData
    {
        public final ThresholdType Type;
        public final double AbsoluteDiff;
        public final double PercentDiff;

        public ThresholdData(final ThresholdType type, final double absoluteDiff, final double percentDiff)
        {
            Type = type;
            AbsoluteDiff = absoluteDiff;
            PercentDiff = percentDiff;
        }

        public boolean hasDiff(double value1, double value2)
        {
            if(value1 == 0 && value2 == 0)
                return false;

            double absDiff = abs(value1 - value2);

            boolean hasAbsDiff = absDiff > AbsoluteDiff;
            boolean hasRelDiff = absDiff / max(value1, value2) > PercentDiff;

            if(Type == ThresholdType.ABSOLUTE_AND_PERCENT)
                return hasAbsDiff && hasRelDiff;

            if(Type == ThresholdType.ABSOLUTE)
                return hasAbsDiff;

            return hasRelDiff;
        }
    }
}
