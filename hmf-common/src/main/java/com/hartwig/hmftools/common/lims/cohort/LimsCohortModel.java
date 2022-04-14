package com.hartwig.hmftools.common.lims.cohort;

import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsCohortModel {

    private static final Logger LOGGER = LogManager.getLogger(LimsCohortModel.class);

    @NotNull
    protected abstract Map<String, LimsCohortConfig> limsCohortMap();

    @Nullable
    public LimsCohortConfig queryCohortData(@Nullable String cohortString, @NotNull String sampleId) {
        String effectiveCohortString = cohortString;
        if (cohortString == null || cohortString.isEmpty()) {
            effectiveCohortString = resolveFromSampleId(sampleId);
            if (effectiveCohortString != null) {
                LOGGER.warn("No cohort string present in LIMS for sample '{}'. Cohort has been set to {}", sampleId, effectiveCohortString);
            }
        }

        LimsCohortConfig cohortConfigData = limsCohortMap().get(effectiveCohortString);
        if (cohortConfigData == null) {
            LOGGER.warn("Could not resolve cohort config for sample '{}' based on LIMS cohort '{}'", sampleId, cohortString);
            return null;
        } else {
            if (cohortConfigData.cohortId().contains("CPCT") || cohortConfigData.cohortId().contains("DRUP")) {
                if (sampleId.startsWith(cohortConfigData.cohortId().substring(0, 4))) {
                    return cohortConfigData;
                } else {
                    throw new IllegalStateException("Sample '" + sampleId + "' does not seem to be part of CPCT/DRUP cohort");
                }
            } else {
                return cohortConfigData;
            }
        }
    }

    @Nullable
    private static String resolveFromSampleId(@NotNull String sampleId) {
        if (sampleId.startsWith("CPCT")) {
            return "CPCT";
        } else if (sampleId.startsWith("DRUP")) {
            return "DRUP";
        } else if (sampleId.startsWith("WIDE")) {
            return "WIDE";
        } else if (sampleId.startsWith("COREDB")) {
            return "COREDB";
        } else if (sampleId.startsWith("CORESC11")) {
            return "CORESC11";
        } else if (sampleId.startsWith("CORELR11")) {
            return "CORELR11";
        } else if (sampleId.startsWith("CORELR02")) {
            return "CORELR02";
        } else if (sampleId.startsWith("CORERI02")) {
            return "CORERI02";
        } else if (sampleId.startsWith("CORE")) {
            return "CORE";
        } else if (sampleId.startsWith("ACTN")) {
            return "ACTIN";
        } else if (sampleId.startsWith("GLOW")) {
            return "GLOW";
        } else if (sampleId.startsWith("OPTC")) {
            return "OPTIC";
        } else if (sampleId.startsWith("SHRP")) {
            return "SHERPA";
        } else if (sampleId.startsWith("GAYA")) {
            return "GENAYA";
        } else if (sampleId.startsWith("OMIC")) {
            return "OMIC";
        } else if (sampleId.startsWith("TARG")) {
            return "TARGTO";
        }else {
            return null;
        }
    }
}