package com.hartwig.hmftools.healthchecker;

import java.util.List;

import com.hartwig.hmftools.common.metrics.WGSMetricQC;
import com.hartwig.hmftools.healthchecker.result.QCValue;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class HealthCheckEvaluation {

    private static final Logger LOGGER = LogManager.getLogger(HealthCheckEvaluation.class);

    private static final String PURPLE_QC_PASS = "PASS";
    private static final String PURPLE_QC_FAIL = "FAIL";
    private static final double MAX_PURPLE_CONTAMINATION = 0.1;
    private static final double MIN_MAPPED_PROPORTION = 0.5;

    private HealthCheckEvaluation() {
    }

    static boolean isPass(@NotNull List<QCValue> qcValues) {
        boolean success = true;
        for (QCValue qcValue : qcValues) {
            if (!succeed(qcValue)) {
                success = false;
            }
        }
        return success;
    }

    private static boolean succeed(@NotNull QCValue qcValue) {
        switch (qcValue.type()) {
            case REF_COVERAGE_10X:
                return checkCoverage(qcValue.value(), "Ref 10x", WGSMetricQC.MIN_REF_10X_COVERAGE);
            case REF_COVERAGE_20X:
                return checkCoverage(qcValue.value(), "Ref 20x", WGSMetricQC.MIN_REF_20X_COVERAGE);
            case TUMOR_COVERAGE_30X:
                return checkCoverage(qcValue.value(), "Tumor 30x", WGSMetricQC.MIN_TUMOR_30X_COVERAGE);
            case TUMOR_COVERAGE_60X:
                return checkCoverage(qcValue.value(), "Tumor 60x", WGSMetricQC.MIN_TUMOR_60X_COVERAGE);
            case PURPLE_QC_STATUS:
                return checkPurpleQCStatus(qcValue.value());
            case PURPLE_CONTAMINATION:
                return checkPurpleContamination(qcValue.value());
            case REF_PROPORTION_MAPPED:
                return checkFlagstatMappingProportion(qcValue.value(), "Ref", MIN_MAPPED_PROPORTION);
            case TUM_PROPORTION_MAPPED:
                return checkFlagstatMappingProportion(qcValue.value(), "Tum", MIN_MAPPED_PROPORTION);
            default: {
                LOGGER.warn("Unrecognized check to evaluate: {}", qcValue.type());
                return false;
            }
        }
    }

    private static boolean checkCoverage(@NotNull String value, @NotNull String name, double minPercentage) {
        double coverage = Double.parseDouble(value);
        if (coverage >= minPercentage) {
            LOGGER.info("PASS - {} coverage percentage of {} is higher than min value {}", name, value, minPercentage);
            return true;
        } else {
            LOGGER.info("FAIL - {} coverage percentage of {} is lower than min value {}", name, value, minPercentage);
            return false;
        }
    }

    private static boolean checkPurpleQCStatus(final String value) {
        if (value.equals(PURPLE_QC_PASS)) {
            LOGGER.info("PASS - Purple QC value is {}", value);
            return true;
        } else if (value.contains(PURPLE_QC_FAIL)) {
            LOGGER.info("FAIL - Purple QC value is {}", value);
            return false;
        } else {
            LOGGER.warn("WARN - Purple QC value is {}", value);
            return true;
        }
    }

    private static boolean checkPurpleContamination(@NotNull String value) {
        double contamination = Double.parseDouble(value);
        if (contamination <= MAX_PURPLE_CONTAMINATION) {
            LOGGER.info("PASS - Contamination of {} is below {}", value, MAX_PURPLE_CONTAMINATION);
            if (contamination > 0) {
                LOGGER.warn("  Contamination found in sample!");
            }
            return true;
        } else {
            LOGGER.info("FAIL - Contamination of {} is above {}", value, MAX_PURPLE_CONTAMINATION);
            return false;
        }
    }

    private static boolean checkFlagstatMappingProportion(@NotNull String value, @NotNull String name, double minProportion) {
        double proportion = Double.parseDouble(value);
        if (proportion >= minProportion) {
            LOGGER.info("PASS - {} proportion of reads mapped {} is higher than min value {}", name, value, minProportion);
            return true;
        } else {
            LOGGER.info("FAIL - {} proportion of reads mapped {} is lower than min value {}", name, value, minProportion);
            return false;
        }
    }
}
