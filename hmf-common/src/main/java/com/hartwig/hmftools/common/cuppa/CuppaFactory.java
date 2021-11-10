package com.hartwig.hmftools.common.cuppa;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CuppaFactory {

    private static final Logger LOGGER = LogManager.getLogger(CuppaFactory.class);

    private static final String SV_TRAIT_CATEGORY = "SV";
    private static final String DNA_COMBINED_DATATYPE = "DNA_COMBINED";

    private CuppaFactory() {
    }

    @NotNull
    public static CuppaData create(@NotNull List<CuppaEntry> cuppaEntries) {
        CuppaEntry best = findMostLikelyPrimaryTumorEntry((cuppaEntries));

        return ImmutableCuppaData.builder()
                .predictedCancerType(best != null ? best.refCancerType() : "Undetermined")
                .bestPredictionLikelihood(best != null ? best.refValue() : 0D)
                .simpleDups32To200B(safeInt(cuppaEntries, "SIMPLE_DUP_32B_200B"))
                .maxComplexSize(safeInt(cuppaEntries, "MAX_COMPLEX_SIZE"))
                .LINECount(safeInt(cuppaEntries, "LINE"))
                .telomericSGLs(safeInt(cuppaEntries, "TELOMERIC_SGL"))
                .build();
    }

    @Nullable
    private static CuppaEntry findMostLikelyPrimaryTumorEntry(@NotNull List<CuppaEntry> cuppaEntries) {
        CuppaEntry best = null;
        for (CuppaEntry entry : cuppaEntries) {
            if (entry.dataType().equals(DNA_COMBINED_DATATYPE)) {
                if (best == null || entry.refValue() > best.refValue()) {
                    best = entry;
                }
            }
        }

        LOGGER.warn("Could not find a single entry of data type '{}'", DNA_COMBINED_DATATYPE);
        return best;
    }

    private static int safeInt(@NotNull List<CuppaEntry> cuppaEntries, @NotNull String dataType) {
        CuppaEntry entry = findEntry(cuppaEntries, dataType);
        if (entry != null) {
            return (int) Math.round(Double.parseDouble(entry.value()));
        } else {
            // -1 is a magic value that can never exist in reality.
            return -1;
        }
    }

    @Nullable
    private static CuppaEntry findEntry(@NotNull List<CuppaEntry> cuppaEntries, @NotNull String dataType) {
        for (CuppaEntry entry : cuppaEntries) {
            if (entry.category().equals(SV_TRAIT_CATEGORY) && entry.dataType().equals(dataType)) {
                return entry;
            }
        }

        LOGGER.warn("Could not find entry with data type '{}'", dataType);
        return null;
    }
}
