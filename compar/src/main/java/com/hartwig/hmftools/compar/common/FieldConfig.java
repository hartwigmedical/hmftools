package com.hartwig.hmftools.compar.common;

import static java.util.Collections.emptyMap;

import java.util.Map;

import com.google.common.collect.Maps;

public class FieldConfig
{
    private final Map<CategoryType, Map<String, ThresholdData>> mFieldThresholds;

    public static final double DEFAULT_DIFF_PERC = 0.1;

    public static final ThresholdData DEFAULT_DECIMAL_THRESHOLD = new ThresholdData(
            ThresholdType.ABSOLUTE_AND_PERCENT, 1, DEFAULT_DIFF_PERC);

    public FieldConfig()
    {
        mFieldThresholds = Maps.newHashMap();
    }

    public boolean isFieldRegistered(final CategoryType category, final String field)
    {
        return mFieldThresholds.containsKey(category) && mFieldThresholds.get(category).containsKey(field);
    }

    public ThresholdData getThreshold(final CategoryType category, final String field)
    {
        return mFieldThresholds.getOrDefault(category, emptyMap()).get(field);
    }

    public boolean hasDifference(final CategoryType category, final String field, double value1, double value2)
    {
        ThresholdData thresholdData = getThreshold(category, field);

        if(thresholdData == null)
            return false;

        return thresholdData.hasDiff(value1, value2);
    }

    public void addFieldThreshold(final CategoryType category, final String field, double absoluteDiff, double percentDiff)
    {
        if(isFieldRegistered(category, field)) // keep any config overrides
            return;

        ThresholdType type = absoluteDiff > 0 && percentDiff > 0 ? ThresholdType.ABSOLUTE_AND_PERCENT :
                (absoluteDiff > 0 ? ThresholdType.ABSOLUTE : ThresholdType.PERCENT);

        mFieldThresholds.putIfAbsent(category, Maps.newHashMap());
        mFieldThresholds.get(category).put(field, new ThresholdData(type, absoluteDiff, percentDiff));
    }
}
