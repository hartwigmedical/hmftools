package com.hartwig.hmftools.healthchecker;

import static com.hartwig.hmftools.common.metrics.WGSMetricQC.REF_COVERAGE_LEVEL_10x;
import static com.hartwig.hmftools.common.metrics.WGSMetricQC.REF_COVERAGE_LEVEL_20x;
import static com.hartwig.hmftools.common.metrics.WGSMetricQC.TUMOR_COVERAGE_LEVEL_30x;
import static com.hartwig.hmftools.common.metrics.WGSMetricQC.TUMOR_COVERAGE_LEVEL_60x;
import static com.hartwig.hmftools.healthchecker.HealthChecksApplication.HC_LOGGER;
import static com.hartwig.hmftools.healthchecker.QCValueType.REF_COVERAGE_10X;
import static com.hartwig.hmftools.healthchecker.QCValueType.REF_COVERAGE_20X;
import static com.hartwig.hmftools.healthchecker.QCValueType.REF_PROPORTION_DUPLICATE;
import static com.hartwig.hmftools.healthchecker.QCValueType.REF_PROPORTION_MAPPED;
import static com.hartwig.hmftools.healthchecker.QCValueType.TUM_COVERAGE_30X;
import static com.hartwig.hmftools.healthchecker.QCValueType.TUM_COVERAGE_60X;
import static com.hartwig.hmftools.healthchecker.QCValueType.TUM_PROPORTION_DUPLICATE;
import static com.hartwig.hmftools.healthchecker.QCValueType.TUM_PROPORTION_MAPPED;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricsSummary;

public class MetricsDataLoader
{
    public static List<QCValue> loadMetricValues(final String sampleId, final String metricsDir, boolean isTumor)
    {
        String summaryFilename = BamMetricsSummary.generateFilename(metricsDir, sampleId);
        String flagStatsFilename = BamFlagStats.generateFilename(metricsDir, sampleId);

        try
        {
            BamMetricsSummary bamMetricsSummary = BamMetricsSummary.read(summaryFilename);
            BamFlagStats bamFlagStats = BamFlagStats.read(flagStatsFilename);

            List<QCValue> qcValues = Lists.newArrayList();

            if(isTumor)
            {
                addCoverageLevel(bamMetricsSummary, TUMOR_COVERAGE_LEVEL_30x, qcValues, TUM_COVERAGE_30X);
                addCoverageLevel(bamMetricsSummary, TUMOR_COVERAGE_LEVEL_60x, qcValues, TUM_COVERAGE_60X);

                qcValues.add(new QCValue(TUM_PROPORTION_MAPPED, String.valueOf(bamFlagStats.mappedProportion())));
                qcValues.add(new QCValue(TUM_PROPORTION_DUPLICATE, String.valueOf(bamFlagStats.duplicateProportion())));
            }
            else
            {
                addCoverageLevel(bamMetricsSummary, REF_COVERAGE_LEVEL_10x, qcValues, REF_COVERAGE_10X);
                addCoverageLevel(bamMetricsSummary, REF_COVERAGE_LEVEL_20x, qcValues, REF_COVERAGE_20X);

                qcValues.add(new QCValue(REF_PROPORTION_MAPPED, String.valueOf(bamFlagStats.mappedProportion())));
                qcValues.add(new QCValue(REF_PROPORTION_DUPLICATE, String.valueOf(bamFlagStats.duplicateProportion())));
            }

            return qcValues;
        }
        catch(IOException e)
        {
            HC_LOGGER.error("failed to read BAM metrics({}): {}", summaryFilename, e.toString());
            return null;
        }
    }

    private static void addCoverageLevel(
            final BamMetricsSummary bamMetricsSummary, int coverageLevel, final List<QCValue> qcValues, final QCValueType qcType)
    {
        Double coveragePercent = bamMetricsSummary.getCoveragePercent(coverageLevel);

        if(coveragePercent != null)
            qcValues.add(new QCValue(qcType, String.valueOf(coveragePercent)));
    }
}
