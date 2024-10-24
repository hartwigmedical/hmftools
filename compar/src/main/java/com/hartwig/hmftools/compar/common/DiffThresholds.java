package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.util.Map;

import com.google.common.collect.Maps;

public class DiffThresholds
{
    private final Map<String, ThresholdData> mFieldThresholds;

    private static final String THRESHOLD_ITEM_DELIM = ":";

    public static final double DEFAULT_DIFF_PERC = 0.1;

    public static final ThresholdData DEFAULT_DECIMAL_THRESHOLD = new ThresholdData(
            ThresholdType.ABSOLUTE_AND_PERCENT, 1, DEFAULT_DIFF_PERC);

    public DiffThresholds()
    {
        mFieldThresholds = Maps.newHashMap();
    }

    public boolean isFieldRegistered(final String field) { return mFieldThresholds.containsKey(field); }

    public boolean hasDifference(final String field, double value1, double value2)
    {
        ThresholdData thresholdData = mFieldThresholds.get(field);

        if(thresholdData == null)
            return false;

        return thresholdData.hasDiff(value1, value2);
    }

    public void addFieldThreshold(final String field, double absoluteDiff, double percentDiff)
    {
        if(mFieldThresholds.containsKey(field)) // keep any config overrides
            return;

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
}
