package com.hartwig.hmftools.patientreporter.algo.orange.cuppa;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.DataTypes;
import com.hartwig.hmftools.common.cuppa.SvDataType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CuppaDataFactory {

    private static final Logger LOGGER = LogManager.getLogger(CuppaDataFactory.class);

    private CuppaDataFactory() {
    }

    @NotNull
    public static CuppaData create(@NotNull List<CuppaDataFile> entries) {
        return ImmutableCuppaData.builder()
                .predictions(extractPredictions(entries))
                .simpleDups32To200B(safeInt(entries, SvDataType.SIMPLE_DUP_32B_200B))
                .maxComplexSize(safeInt(entries, SvDataType.MAX_COMPLEX_SIZE))
                .LINECount(safeInt(entries, SvDataType.LINE))
                .telomericSGLs(safeInt(entries, SvDataType.TELOMERIC_SGL))
                .build();
    }

    @NotNull
    private static List<CuppaPrediction> extractPredictions(@NotNull List<CuppaDataFile> entries) {
        String bestCombinedType = determineBestCombinedDataType(entries);
        if (bestCombinedType == null) {
            LOGGER.warn("Could not find a valid combined data type amongst cuppa entries");
            return Lists.newArrayList();
        }

        List<CuppaPrediction> predictions = Lists.newArrayList();

        for (CuppaDataFile entry : entries) {
            if (entry.DataType.equals(bestCombinedType)) {
                for (Map.Entry<String, Double> cancerTypeEntry : entry.CancerTypeValues.entrySet()) {
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

        for (CuppaDataFile entry : entries) {
            if (entry.Category == CategoryType.COMBINED) {
                switch (entry.DataType) {
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

        if (hasOverallCombinedType) {
            return DataTypes.DATA_TYPE_COMBINED;
        } else if (hasDnaCombinedType) {
            return DataTypes.DATA_TYPE_DNA_COMBINED;
        } else if (hasRnaCombinedType) {
            return DataTypes.DATA_TYPE_RNA_COMBINED;
        }

        return null;
    }

    private static int safeInt(@NotNull List<CuppaDataFile> entries, @NotNull SvDataType svDataType) {
        CuppaDataFile entry = findSvEntry(entries, svDataType);
        if (entry != null) {
            return (int) Math.round(Double.parseDouble(entry.Value));
        } else {
            // -1 is a magic value that can never exist in reality.
            return -1;
        }
    }

    @Nullable
    private static CuppaDataFile findSvEntry(@NotNull List<CuppaDataFile> entries, @NotNull SvDataType dataType) {
        for (CuppaDataFile entry : entries) {
            if (entry.Category == CategoryType.SV && entry.DataType.equals(dataType.toString())) {
                return entry;
            }
        }

        LOGGER.warn("Could not find entry with data type '{}'", dataType);
        return null;
    }
}
