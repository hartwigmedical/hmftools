package com.hartwig.hmftools.healthchecker;

import static com.hartwig.hmftools.common.metrics.WGSMetricQC.MIN_MAPPED_PROPORTION;
import static com.hartwig.hmftools.common.metrics.WGSMetricQC.hasSufficientMappedProportion;
import static com.hartwig.hmftools.healthchecker.HealthChecksApplication.HC_LOGGER;

import java.util.List;

import com.hartwig.hmftools.common.metrics.WGSMetricQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;

final class HealthCheckEvaluation
{
    private static final String PURPLE_QC_PASS = "PASS";
    private static final String PURPLE_QC_FAIL = "FAIL";

    static boolean isPass(final List<QCValue> qcValues)
    {
        for(QCValue qcValue : qcValues)
        {
            if(!succeed(qcValue))
                return false;
        }

        return true;
    }

    private static boolean succeed(final QCValue qcValue)
    {
        switch(qcValue.Type)
        {
            case REF_COVERAGE_10X:
                return checkCoverage(qcValue.Value, "Ref 10x", WGSMetricQC.MIN_REF_10X_COVERAGE);
            case REF_COVERAGE_20X:
                return checkCoverage(qcValue.Value, "Ref 20x", WGSMetricQC.MIN_REF_20X_COVERAGE);
            case TUM_COVERAGE_30X:
                return checkCoverage(qcValue.Value, "Tum 30x", WGSMetricQC.MIN_TUMOR_30X_COVERAGE);
            case TUM_COVERAGE_60X:
                return checkCoverage(qcValue.Value, "Tum 60x", WGSMetricQC.MIN_TUMOR_60X_COVERAGE);
            case PURPLE_QC_STATUS:
                return checkPurpleQCStatus(qcValue.Value);
            case PURPLE_CONTAMINATION:
                return checkPurpleContamination(qcValue.Value);
            case REF_PROPORTION_MAPPED:
                return checkFlagstatMappingProportion(qcValue.Value, "Ref");
            case TUM_PROPORTION_MAPPED:
                return checkFlagstatMappingProportion(qcValue.Value, "Tum");
            case REF_PROPORTION_DUPLICATE:
            case TUM_PROPORTION_DUPLICATE:
                // No QC on duplicate rate (only reporting the value)
                return true;
            default:
            {
                HC_LOGGER.warn("Unrecognized check to evaluate: {}", qcValue.Type);
                return false;
            }
        }
    }

    private static boolean checkCoverage(final String value, final String name, double minPercentage)
    {
        double coverage = Double.parseDouble(value);
        if(coverage >= minPercentage)
        {
            HC_LOGGER.info("QC PASS - {} coverage of {} is higher than min value {}", name, value, minPercentage);
            return true;
        }
        else
        {
            HC_LOGGER.info("QC FAIL - {} coverage of {} is lower than min value {}", name, value, minPercentage);
            return false;
        }
    }

    private static boolean checkPurpleQCStatus(final String value)
    {
        if(value.equals(PURPLE_QC_PASS))
        {
            HC_LOGGER.info("QC PASS - Purple QC value is {}", value);
            return true;
        }
        else if(value.contains(PURPLE_QC_FAIL))
        {
            HC_LOGGER.info("QC FAIL - Purple QC value is {}", value);
            return false;
        }
        else
        {
            HC_LOGGER.warn("QC WARN - Purple QC value is {}", value);
            return true;
        }
    }

    private static boolean checkPurpleContamination(final String value)
    {
        double contamination = Double.parseDouble(value);
        if(contamination <= PurpleQCStatus.MAX_CONTAMINATION)
        {
            HC_LOGGER.info("QC PASS - Contamination of {} is lower than {}", value, PurpleQCStatus.MAX_CONTAMINATION);
            if(contamination > 0)
            {
                HC_LOGGER.warn("  But contamination is higher than 0!");
            }
            return true;
        }
        else
        {
            HC_LOGGER.info("QC FAIL - Contamination of {} is higher than {}", value, PurpleQCStatus.MAX_CONTAMINATION);
            return false;
        }
    }

    private static boolean checkFlagstatMappingProportion(final String value, final String name)
    {
        double proportion = Double.parseDouble(value);

        if(hasSufficientMappedProportion(proportion))
        {
            HC_LOGGER.info("QC PASS - {} mapping percentage {} is higher than min value {}", name, value, MIN_MAPPED_PROPORTION);
            return true;
        }
        else
        {
            HC_LOGGER.info("QC FAIL - {} mapping percentage {} is lower than min value {}", name, value, MIN_MAPPED_PROPORTION);
            return false;
        }
    }
}
