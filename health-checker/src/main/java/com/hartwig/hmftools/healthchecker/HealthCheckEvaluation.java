package com.hartwig.hmftools.healthchecker;

import java.util.List;

import com.hartwig.hmftools.healthchecker.result.QCValue;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class HealthCheckEvaluation {

    private static final Logger LOGGER = LogManager.getLogger(HealthCheckEvaluation.class);

    private static final double MIN_REF_10X_COVERAGE = 0.9;
    private static final double MIN_REF_20X_COVERAGE = 0.7;
    private static final double MIN_TUMOR_30X_COVERAGE = 0.8;
    private static final double MIN_TUMOR_60X_COVERAGE = 0.65;

    // Complete sample swap  =>  <0.4
    // 10%+ contamination in tumor =>  < 0.486
    // Potential Tumor contamination in Ref  => >0.5
    private static final double MIN_AMBER_MEAN_BAF = 0.48;
    private static final double MAX_AMBER_MEAN_BAF = 0.6;

    private static final double MAX_AMBER_CONTAMINATION = 0.1;

    private static final String PURPLE_QC_PASS = "PASS";

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
                return checkCoverage(qcValue.value(), "Ref 10x", MIN_REF_10X_COVERAGE);
            case REF_COVERAGE_20X:
                return checkCoverage(qcValue.value(), "Ref 20x", MIN_REF_20X_COVERAGE);
            case TUMOR_COVERAGE_30X:
                return checkCoverage(qcValue.value(), "Tumor 30x", MIN_TUMOR_30X_COVERAGE);
            case TUMOR_COVERAGE_60X:
                return checkCoverage(qcValue.value(), "Tumor 60x", MIN_TUMOR_60X_COVERAGE);
            case AMBER_MEAN_BAF:
                return checkAmberMeanBAF(qcValue.value());
            case AMBER_CONTAMINATION:
                return checkAmberContamination(qcValue.value());
            case PURPLE_QC_STATUS:
                return checkPurpleQCStatus(qcValue.value());
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

    private static boolean checkAmberMeanBAF(@NotNull String value) {
        double meanBaf = Double.parseDouble(value);

        if (meanBaf >= MIN_AMBER_MEAN_BAF && meanBaf <= MAX_AMBER_MEAN_BAF) {
            LOGGER.info("PASS - Mean BAF of {} is between {} and {}", value, MIN_AMBER_MEAN_BAF, MAX_AMBER_MEAN_BAF);
            return true;
        } else if (meanBaf < MIN_AMBER_MEAN_BAF) {
            LOGGER.info("FAIL - Mean BAF of {} is below {}", value, MIN_AMBER_MEAN_BAF);
            return false;
        } else {
            LOGGER.info("FAIL - Mean BAF of {} is above {}", value, MAX_AMBER_MEAN_BAF);
            return false;
        }
    }

    private static boolean checkAmberContamination(@NotNull String value) {
        double contamination = Double.parseDouble(value);

        if (contamination <= MAX_AMBER_CONTAMINATION) {
            LOGGER.info("PASS - Contamination of {} is below {}", value, MAX_AMBER_CONTAMINATION);
            if (contamination > 0) {
                LOGGER.warn("  Contamination found in sample!");
            }
            return true;
        } else {
            LOGGER.info("FAIL - Contamination of {} is above {}", value, MAX_AMBER_CONTAMINATION);
            return false;
        }
    }

    private static boolean checkPurpleQCStatus(final String value) {
        if (value.equals(PURPLE_QC_PASS)) {
            LOGGER.info("PASS - Purple QC value is {}", value);
            return true;
        } else {
            LOGGER.info("FAIL - Purple QC value is {}", value);
            return false;
        }
    }
}
