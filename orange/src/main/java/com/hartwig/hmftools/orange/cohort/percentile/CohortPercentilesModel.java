package com.hartwig.hmftools.orange.cohort.percentile;

import java.util.Collection;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.orange.cohort.datamodel.Evaluation;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableEvaluation;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;
import com.hartwig.hmftools.orange.cohort.mapping.CohortConstants;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapper;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CohortPercentilesModel {

    private static final Logger LOGGER = LogManager.getLogger(CohortPercentilesModel.class);

    @NotNull
    private final CohortMapper mapper;
    @NotNull
    private final Multimap<PercentileType, CohortPercentiles> percentileMap;

    public CohortPercentilesModel(@NotNull final CohortMapper mapper,
            @NotNull final Multimap<PercentileType, CohortPercentiles> percentileMap) {
        this.mapper = mapper;
        this.percentileMap = percentileMap;
    }

    @Nullable
    public Evaluation percentile(@NotNull Observation observation) {
        Collection<CohortPercentiles> percentiles = percentileMap.get(observation.type());
        if (percentiles == null) {
            LOGGER.warn("Percentile evaluation requested for a type that is not supported by model: {}", observation.type());
            return null;
        }

        double panCancerPercentile = determinePercentile(percentiles, CohortConstants.COHORT_PAN_CANCER, observation.value());
        Double cancerTypePercentile = null;
        String cancerType = mapper.cancerTypeForSample(observation.sample());
        if (cancerType != null) {
            cancerTypePercentile = determinePercentile(percentiles, cancerType, observation.value());
        }

        return ImmutableEvaluation.builder()
                .cancerType(cancerType)
                .panCancerPercentile(panCancerPercentile)
                .cancerTypePercentile(cancerTypePercentile)
                .build();
    }

    private double determinePercentile(@NotNull Collection<CohortPercentiles> percentiles, @NotNull String cohort, double value) {
        CohortPercentiles percentile = find(percentiles, cohort);

        if (percentile == null) {
            throw new IllegalStateException("Could not find cohort '" + cohort + "' in list of percentiles!");
        }

        int index = 0;
        double refValue = percentile.values().get(index);
        while (value > refValue && index < percentile.values().size()) {
            refValue = percentile.values().get(index);
            index++;
        }

        return (double) index / percentile.values().size();
    }

    @Nullable
    private static CohortPercentiles find(@NotNull Collection<CohortPercentiles> percentiles, @NotNull String cohort) {
        for (CohortPercentiles percentile : percentiles) {
            if (percentile.cancerType().equals(cohort)) {
                return percentile;
            }
        }
        return null;
    }
}
