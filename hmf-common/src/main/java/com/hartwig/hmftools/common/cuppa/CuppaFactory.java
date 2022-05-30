package com.hartwig.hmftools.common.cuppa;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CuppaFactory {

    private static final Logger LOGGER = LogManager.getLogger(CuppaFactory.class);

    private static final String SV_TRAIT_CATEGORY = "SV";
    static final String COMBINED_CATEGORY = "COMBINED";

    static final String DNA_COMBINED_DATATYPE = "DNA_COMBINED";
    static final String RNA_COMBINED_DATATYPE = "RNA_COMBINED";
    static final String OVERALL_COMBINED_DATATYPE = "COMBINED";

    private CuppaFactory() {
    }

    @NotNull
    public static CuppaData create(@NotNull List<CuppaEntry> entries) {
        return ImmutableCuppaData.builder()
                .predictions(extractPredictions(entries))
                .simpleDups32To200B(safeInt(entries, "SIMPLE_DUP_32B_200B"))
                .maxComplexSize(safeInt(entries, "MAX_COMPLEX_SIZE"))
                .LINECount(safeInt(entries, "LINE"))
                .telomericSGLs(safeInt(entries, "TELOMERIC_SGL"))
                .build();
    }

    @NotNull
    private static List<CuppaPrediction> extractPredictions(@NotNull List<CuppaEntry> entries) {
        String bestCombinedType = determineBestCombinedDataType(entries);
        if (bestCombinedType == null) {
            LOGGER.warn("Could not find a valid combined data type amongst cuppa entries");
            return Lists.newArrayList();
        }

        List<CuppaPrediction> predictions = Lists.newArrayList();

        for (CuppaEntry entry : entries) {
            if (entry.dataType().equals(bestCombinedType)) {
                predictions.add(ImmutableCuppaPrediction.builder()
                        .cancerType(entry.refCancerType())
                        .likelihood(entry.refValue())
                        .build());
            }
        }

        predictions.sort(new CuppaPredictionComparator());

        return predictions;
    }

    @Nullable
    private static String determineBestCombinedDataType(@NotNull List<CuppaEntry> entries) {
        boolean hasDnaCombinedType = false;
        boolean hasRnaCombinedType = false;
        boolean hasOverallCombinedType = false;

        for (CuppaEntry entry : entries) {
            if (entry.category().equals(COMBINED_CATEGORY)) {
                if (entry.dataType().equals(OVERALL_COMBINED_DATATYPE)) {
                    hasOverallCombinedType = true;
                } else if (entry.dataType().equals(DNA_COMBINED_DATATYPE)) {
                    hasDnaCombinedType = true;
                } else if (entry.dataType().equals(RNA_COMBINED_DATATYPE)) {
                    hasRnaCombinedType = true;
                } else {
                    LOGGER.warn("Unrecognized combined data type: {}", entry.dataType());
                }
            }
        }

        if (hasOverallCombinedType) {
            return OVERALL_COMBINED_DATATYPE;
        } else if (hasDnaCombinedType) {
            return DNA_COMBINED_DATATYPE;
        } else if (hasRnaCombinedType) {
            return RNA_COMBINED_DATATYPE;
        }

        return null;
    }

    private static int safeInt(@NotNull List<CuppaEntry> entries, @NotNull String dataType) {
        CuppaEntry entry = findSvEntry(entries, dataType);
        if (entry != null) {
            return (int) Math.round(Double.parseDouble(entry.value()));
        } else {
            // -1 is a magic value that can never exist in reality.
            return -1;
        }
    }

    @Nullable
    private static CuppaEntry findSvEntry(@NotNull List<CuppaEntry> entries, @NotNull String dataType) {
        for (CuppaEntry entry : entries) {
            if (entry.category().equals(SV_TRAIT_CATEGORY) && entry.dataType().equals(dataType)) {
                return entry;
            }
        }

        LOGGER.warn("Could not find entry with data type '{}'", dataType);
        return null;
    }
}
