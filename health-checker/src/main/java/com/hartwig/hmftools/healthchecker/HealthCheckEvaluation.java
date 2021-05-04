package com.hartwig.hmftools.healthchecker;

import java.util.List;

import com.hartwig.hmftools.common.flagstat.FlagstatQC;
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
            case TUM_COVERAGE_30X:
                return checkCoverage(qcValue.value(), "Tum 30x", WGSMetricQC.MIN_TUMOR_30X_COVERAGE);
            case TUM_COVERAGE_60X:
                return checkCoverage(qcValue.value(), "Tum 60x", WGSMetricQC.MIN_TUMOR_60X_COVERAGE);
            case PURPLE_QC_STATUS:
                return checkPurpleQCStatus(qcValue.value());
            case PURPLE_CONTAMINATION:
                return checkPurpleContamination(qcValue.value());
            case REF_PROPORTION_MAPPED:
                return checkFlagstatMappingProportion(qcValue.value(), "Ref");
            case TUM_PROPORTION_MAPPED:
                return checkFlagstatMappingProportion(qcValue.value(), "Tum");
            case REF_PROPORTION_DUPLICATE:
            case TUM_PROPORTION_DUPLICATE:
                // No QC on duplicate rate (only reporting the value)
                return true;
            default: {
                LOGGER.warn("Unrecognized check to evaluate: {}", qcValue.type());
                return false;
            }
        }
    }

    private static boolean checkCoverage(@NotNull String value, @NotNull String name, double minPercentage) {
        double coverage = Double.parseDouble(value);
        if (coverage >= minPercentage) {
            LOGGER.info("QC PASS - {} coverage of {} is higher than min value {}", name, value, minPercentage);
            return true;
        } else {
            LOGGER.info("QC FAIL - {} coverage of {} is lower than min value {}", name, value, minPercentage);
            return false;
        }
    }

    private static boolean checkPurpleQCStatus(@NotNull final String value) {
        if (value.equals(PURPLE_QC_PASS)) {
            LOGGER.info("QC PASS - Purple QC value is {}", value);
            return true;
        } else if (value.contains(PURPLE_QC_FAIL)) {
            LOGGER.info("QC FAIL - Purple QC value is {}", value);
            return false;
        } else {
            LOGGER.warn("QC WARN - Purple QC value is {}", value);
            return true;
        }
    }

    private static boolean checkPurpleContamination(@NotNull String value) {
        double contamination = Double.parseDouble(value);
        if (contamination <= MAX_PURPLE_CONTAMINATION) {
            LOGGER.info("QC PASS - Contamination of {} is lower than {}", value, MAX_PURPLE_CONTAMINATION);
            if (contamination > 0) {
                LOGGER.warn("  But contamination is higher than 0!");
            }
            return true;
        } else {
            LOGGER.info("QC FAIL - Contamination of {} is higher than {}", value, MAX_PURPLE_CONTAMINATION);
            return false;
        }
    }

    private static boolean checkFlagstatMappingProportion(@NotNull String value, @NotNull String name) {
        double proportion = Double.parseDouble(value);
        if (FlagstatQC.pass(proportion)) {
            LOGGER.info("QC PASS - {} mapping percentage {} is higher than min value {}", name, value, FlagstatQC.MIN_MAPPED_PROPORTION);
            return true;
        } else {
            LOGGER.info("QC FAIL - {} mapping percentage {} is lower than min value {}", name, value, FlagstatQC.MIN_MAPPED_PROPORTION);
            return false;
        }
    }
}
