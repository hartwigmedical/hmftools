package com.hartwig.hmftools.qsee.prep.category;

import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_10;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_100;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_20;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_250;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_30;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_60;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.DUAL_STRAND_READS;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.DUPLICATE_READS;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.LOW_BASE_QUAL;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.LOW_MAP_QUAL;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.MAPPED_PROPORTION;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.MEAN_COVERAGE;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.putFeature;

import java.io.IOException;
import java.util.EnumMap;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CategoryPrepTask;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;
import com.hartwig.hmftools.qsee.prep.category.table.BamMetricsData;
import com.hartwig.hmftools.qsee.status.QcStatus;
import com.hartwig.hmftools.qsee.status.QcStatusType;
import com.hartwig.hmftools.qsee.status.QcThreshold;
import com.hartwig.hmftools.qsee.status.ThresholdRegistry;

import org.jetbrains.annotations.NotNull;

public class SummaryTableBamMetricsPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    public SummaryTableBamMetricsPrep(CommonPrepConfig config) { mConfig = config; }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamMetricSummary bamMetricSummary,
            SampleType sampleType, ThresholdRegistry qcThresholds)
    {
        if(bamMetricSummary == null)
            return;

        putFeature(featuresMap, MEAN_COVERAGE, bamMetricSummary.meanCoverage(),
                qcThresholds.getThreshold(sampleType, MEAN_COVERAGE, QcStatusType.WARN));

        putFeature(featuresMap, LOW_MAP_QUAL, bamMetricSummary.lowMapQualPercent(),
                qcThresholds.getThreshold(sampleType, LOW_MAP_QUAL, QcStatusType.WARN));

        putFeature(featuresMap, LOW_BASE_QUAL, bamMetricSummary.lowBaseQualPercent(),
                qcThresholds.getThreshold(sampleType, LOW_BASE_QUAL, QcStatusType.WARN));

        putFeature(featuresMap, DUPLICATE_READS, (double) bamMetricSummary.duplicateReads() / bamMetricSummary.totalReads(),
                qcThresholds.getThreshold(sampleType, DUPLICATE_READS, QcStatusType.WARN));

        putFeature(featuresMap, DUAL_STRAND_READS, (double) bamMetricSummary.dualStrandReads() / bamMetricSummary.totalReads(),
                qcThresholds.getThreshold(sampleType, DUAL_STRAND_READS, QcStatusType.WARN));
    }

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamMetricCoverage bamMetricCoverage,
            SampleType sampleType, ThresholdRegistry qcThresholds)
    {
        if(bamMetricCoverage == null)
            return;

        List<SummaryTableFeature> coverageAboveXFeatures = List.of(
                COVERAGE_ABOVE_10, COVERAGE_ABOVE_20, COVERAGE_ABOVE_30, COVERAGE_ABOVE_60, COVERAGE_ABOVE_100, COVERAGE_ABOVE_250);

        for(SummaryTableFeature coverageAboveXFeature : coverageAboveXFeatures)
        {
            int coverageThreshold = switch(coverageAboveXFeature)
            {
                case COVERAGE_ABOVE_10 -> 10;
                case COVERAGE_ABOVE_20 -> 20;
                case COVERAGE_ABOVE_30 -> 30;
                case COVERAGE_ABOVE_60 -> 60;
                case COVERAGE_ABOVE_100 -> 100;
                case COVERAGE_ABOVE_250 -> 250;
                default -> throw new IllegalStateException("Unexpected min coverage feature: " + coverageAboveXFeature);
            };

            List<ValueFrequency> coverageBaseCounts = bamMetricCoverage.Coverage;

            long totalBases = coverageBaseCounts.stream().mapToLong(x -> x.Count).sum();

            long basesAboveCoverage = coverageBaseCounts.stream()
                    .filter(x -> x.Value >= coverageThreshold)
                    .mapToLong(x -> x.Count)
                    .sum();

            double propBasesAboveCoverage = (double) basesAboveCoverage / totalBases;

            QcThreshold qcThreshold = qcThresholds.getThreshold(sampleType, coverageAboveXFeature, QcStatusType.WARN);
            QcStatus qcStatus = qcThreshold.getQcStatus(propBasesAboveCoverage);

            putFeature(featuresMap, coverageAboveXFeature, propBasesAboveCoverage, qcStatus);
        }
    }

    private static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamFlagStats bamFlagStats,
            SampleType sampleType, ThresholdRegistry qcThresholds)
    {
        if(bamFlagStats == null)
            return;

        QcStatus qcStatus = qcThresholds
                .getThreshold(sampleType, MAPPED_PROPORTION, QcStatusType.FAIL)
                .getQcStatus(bamFlagStats.mappedProportion());

        putFeature(featuresMap, MAPPED_PROPORTION, bamFlagStats.mappedProportion(), qcStatus);
    }

    @VisibleForTesting
    static List<Feature> createFeatures(BamMetricsData bamMetricsData, SampleType sampleType, ThresholdRegistry qcThresholds)
    {
        EnumMap<SummaryTableFeature, Feature> featuresMap = new EnumMap<>(SummaryTableFeature.class);

        putFeatures(featuresMap, bamMetricsData.bamMetricSummary(), sampleType, qcThresholds);
        putFeatures(featuresMap, bamMetricsData.bamMetricCoverage(), sampleType, qcThresholds);
        putFeatures(featuresMap, bamMetricsData.bamFlagStats(), sampleType, qcThresholds);

        return featuresMap.values().stream().toList();
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        BamMetricsData bamMetricsData = new BamMetricsData(mConfig, sampleId, sampleType);

        if(!bamMetricsData.missingInputPaths().isEmpty())
        {
            CategoryPrepTask.missingInputFilesError(mConfig.AllowMissingInput, this, sampleType, bamMetricsData.formMissingInputsString());
        }

        List<Feature> features = createFeatures(bamMetricsData, sampleType, mConfig.QcThresholds);

        return features;
    }

}
