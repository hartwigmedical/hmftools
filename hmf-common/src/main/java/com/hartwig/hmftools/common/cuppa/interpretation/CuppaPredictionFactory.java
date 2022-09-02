package com.hartwig.hmftools.common.cuppa.interpretation;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.DataTypes;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CuppaPredictionFactory
{
    private static final Logger LOGGER = LogManager.getLogger(CuppaPredictionFactory.class);

    private CuppaPredictionFactory()
    {
    }

    @NotNull
    public static List<CuppaPrediction> create(@NotNull List<CuppaDataFile> entries)
    {
        String bestCombinedType = determineBestCombinedDataType(entries);
        if (bestCombinedType == null)
        {
            LOGGER.warn("Could not find a valid combined data type amongst cuppa entries");
            return Lists.newArrayList();
        }

        List<CuppaPrediction> predictions = Lists.newArrayList();

        for (CuppaDataFile entry : entries)
        {
            if (entry.DataType.equals(bestCombinedType))
            {
                for (Map.Entry<String, Double> cancerTypeEntry : entry.CancerTypeValues.entrySet())
                {
                    predictions.add(ImmutableCuppaPrediction.builder()
                            .cancerType(cancerTypeEntry.getKey())
                            .likelihood(cancerTypeEntry.getValue())
                            .build());
                }
            }
        }

        predictions.sort(new CuppaPredictionComparator());

        return predictions;
    }

    @Nullable
    private static String determineBestCombinedDataType(@NotNull List<CuppaDataFile> entries) {
        boolean hasDnaCombinedType = false;
        boolean hasRnaCombinedType = false;
        boolean hasOverallCombinedType = false;

        for (CuppaDataFile entry : entries)
        {
            if (entry.Category == CategoryType.COMBINED)
            {
                switch (entry.DataType)
                {
                    case DataTypes.DATA_TYPE_COMBINED:
                        hasOverallCombinedType = true;
                        break;
                    case DataTypes.DATA_TYPE_DNA_COMBINED:
                        hasDnaCombinedType = true;
                        break;
                    case DataTypes.DATA_TYPE_RNA_COMBINED:
                        hasRnaCombinedType = true;
                        break;
                    default:
                        LOGGER.warn("Unrecognized combined data type: {}", entry.DataType);
                        break;
                }
            }
        }

        if (hasOverallCombinedType)
        {
            return DataTypes.DATA_TYPE_COMBINED;
        }
        else if (hasDnaCombinedType)
        {
            return DataTypes.DATA_TYPE_DNA_COMBINED;
        }
        else if (hasRnaCombinedType)
        {
            return DataTypes.DATA_TYPE_RNA_COMBINED;
        }

        return null;
    }
}
