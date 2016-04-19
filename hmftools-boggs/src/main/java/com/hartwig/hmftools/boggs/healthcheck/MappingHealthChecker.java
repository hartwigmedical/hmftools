package com.hartwig.hmftools.boggs.healthcheck;

import com.hartwig.hmftools.boggs.PatientData;
import com.hartwig.hmftools.boggs.SampleData;
import com.hartwig.hmftools.boggs.flagstatreader.FlagStatData;
import com.hartwig.hmftools.boggs.flagstatreader.FlagStats;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class MappingHealthChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(MappingHealthChecker.class);

    private static final double MIN_MAPPED_PERCENTAGE = 0.992;
    private static final double MIN_PROPERLY_PAIRED_PERCENTAGE = 0.99;
    private static final double MAX_SINGLETONS = 0.005;
    private static final double MAX_MATE_MAPPED_TO_DIFFERENT_CHR = 0.0001;

    public boolean isHealthy(@NotNull PatientData patient) {
        checkSample(patient.refSample());
        checkSample(patient.tumorSample());
        return true;
    }

    private void checkSample(@NotNull SampleData sample) {
        LOGGER.info("Checking mapping health for " + sample.externalID());
        for (FlagStatData flagstatData : sample.rawMappingFlagstats()) {
            LOGGER.info(" Verifying " + flagstatData.path());
            FlagStats passed = flagstatData.qcPassedReads();
            double mappedPercentage = passed.mapped() / (double) passed.total();
            double properlyPairedPercentage = passed.properlyPaired() / (double) passed.total();
            double singletonPercentage = passed.singletons() / (double) passed.total();
            double mateMappedToDifferentChrPercentage = passed.mateMappedToDifferentChr() / (double) passed.total();

            if (mappedPercentage < MIN_MAPPED_PERCENTAGE) {
                LOGGER.info("  WARN: Low mapped percentage: " +
                        toPercentage(mappedPercentage));
            } else {
                LOGGER.info("  OK: Acceptable mapped percentage: " +
                        toPercentage(mappedPercentage));
            }

            if (properlyPairedPercentage < MIN_PROPERLY_PAIRED_PERCENTAGE) {
                LOGGER.info("  WARN: Low properly paired percentage: " +
                        toPercentage(properlyPairedPercentage));
            } else {
                LOGGER.info("  OK: Acceptable properly paired percentage: " +
                        toPercentage(properlyPairedPercentage));
            }

            if (singletonPercentage > MAX_SINGLETONS) {
                LOGGER.info("  WARN: High singleton percentage: " +
                        toPercentage(singletonPercentage));
            } else {
                LOGGER.info("  OK: Acceptable singleton percentage: " +
                        toPercentage(singletonPercentage));
            }

            if (mateMappedToDifferentChrPercentage > MAX_MATE_MAPPED_TO_DIFFERENT_CHR) {
                LOGGER.info("  WARN: High mate mapped to different chr percentage: " +
                        toPercentage(mateMappedToDifferentChrPercentage));
            } else {
                LOGGER.info("  OK: Acceptable mate mapped to different chr percentage: " +
                        toPercentage(mateMappedToDifferentChrPercentage));
            }
        }
    }

    @NotNull
    private static String toPercentage(double percentage){
        return (Math.round(percentage * 10000L) / 100D) + "%";
    }
}
