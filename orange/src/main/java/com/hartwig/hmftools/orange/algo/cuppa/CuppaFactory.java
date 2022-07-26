package com.hartwig.hmftools.orange.algo.cuppa;

import static com.hartwig.hmftools.common.cuppa.CategoryType.COMBINED;
import static com.hartwig.hmftools.common.cuppa.CategoryType.SV;
import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_COMBINED;
import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_DNA_COMBINED;
import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_RNA_COMBINED;
import static com.hartwig.hmftools.common.cuppa.SvDataType.LINE;
import static com.hartwig.hmftools.common.cuppa.SvDataType.MAX_COMPLEX_SIZE;
import static com.hartwig.hmftools.common.cuppa.SvDataType.SIMPLE_DUP_32B_200B;
import static com.hartwig.hmftools.common.cuppa.SvDataType.TELOMERIC_SGL;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CuppaEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CuppaFactory {

    private static final Logger LOGGER = LogManager.getLogger(CuppaFactory.class);

    private CuppaFactory() {
    }

    @NotNull
    public static CuppaData create(@NotNull List<CuppaEntry> entries) {
        return ImmutableCuppaData.builder()
                .predictions(extractPredictions(entries))
                .simpleDups32To200B(safeInt(entries, SIMPLE_DUP_32B_200B.toString()))
                .maxComplexSize(safeInt(entries, MAX_COMPLEX_SIZE.toString()))
                .LINECount(safeInt(entries, LINE.toString()))
                .telomericSGLs(safeInt(entries, TELOMERIC_SGL.toString()))
                .build();
    }

    private static List<CuppaPrediction> extractPredictions(final List<CuppaEntry> entries) {
        String bestCombinedType = determineBestCombinedDataType(entries);
        if (bestCombinedType == null) {
            LOGGER.warn("Could not find a valid combined data type amongst cuppa entries");
            return Lists.newArrayList();
        }

        List<CuppaPrediction> predictions = Lists.newArrayList();

        for (CuppaEntry entry : entries) {
            if (entry.dataType().equals(bestCombinedType)) {
                predictions.add(ImmutableCuppaPrediction.builder().cancerType(entry.refCancerType()).likelihood(entry.refValue()).build());
            }
        }

        predictions.sort(new CuppaPredictionComparator());

        return predictions;
    }

    private static String determineBestCombinedDataType(final List<CuppaEntry> entries) {
        boolean hasDnaCombinedType = false;
        boolean hasRnaCombinedType = false;
        boolean hasOverallCombinedType = false;

        for (CuppaEntry entry : entries) {
            if (entry.category().equals(COMBINED)) {
                if (entry.dataType().equals(DATA_TYPE_COMBINED)) {
                    hasOverallCombinedType = true;
                } else if (entry.dataType().equals(DATA_TYPE_DNA_COMBINED)) {
                    hasDnaCombinedType = true;
                } else if (entry.dataType().equals(DATA_TYPE_RNA_COMBINED)) {
                    hasRnaCombinedType = true;
                } else {
                    LOGGER.warn("Unrecognized combined data type: {}", entry.dataType());
                }
            }
        }

        if (hasOverallCombinedType) {
            return DATA_TYPE_COMBINED;
        } else if (hasDnaCombinedType) {
            return DATA_TYPE_DNA_COMBINED;
        } else if (hasRnaCombinedType) {
            return DATA_TYPE_RNA_COMBINED;
        }

        return null;
    }

    private static int safeInt(final List<CuppaEntry> entries, final String dataType) {
        CuppaEntry entry = findSvEntry(entries, dataType);
        if (entry != null) {
            return (int) Math.round(Double.parseDouble(entry.value()));
        } else {
            // -1 is a magic value that can never exist in reality.
            return -1;
        }
    }

    @Nullable
    private static CuppaEntry findSvEntry(final List<CuppaEntry> entries, final String dataType) {
        for (CuppaEntry entry : entries) {
            if (entry.category() == SV && entry.dataType().equals(dataType)) {
                return entry;
            }
        }

        LOGGER.warn("Could not find entry with data type '{}'", dataType);
        return null;
    }
}
