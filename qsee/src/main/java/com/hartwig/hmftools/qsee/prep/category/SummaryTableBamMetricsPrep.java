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
import com.hartwig.hmftools.qsee.status.QcThreshold;
import com.hartwig.hmftools.qsee.status.QcThresholdRegistry;
import com.hartwig.hmftools.qsee.status.ThresholdGroup;

import org.jetbrains.annotations.NotNull;

public class SummaryTableBamMetricsPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    public SummaryTableBamMetricsPrep(CommonPrepConfig config) { mConfig = config; }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamMetricSummary bamMetricSummary, SampleType sampleType)
    {
        if(bamMetricSummary == null)
            return;

        putFeature(featuresMap, MEAN_COVERAGE, bamMetricSummary.meanCoverage(),
                QcThresholdRegistry.getThreshold(sampleType, MEAN_COVERAGE));

        putFeature(featuresMap, LOW_MAP_QUAL, bamMetricSummary.lowMapQualPercent(),
                QcThresholdRegistry.getCommonThreshold(LOW_MAP_QUAL));

        putFeature(featuresMap, LOW_BASE_QUAL, bamMetricSummary.lowBaseQualPercent(),
                QcThresholdRegistry.getCommonThreshold(LOW_BASE_QUAL));

        putFeature(featuresMap, DUPLICATE_READS, (double) bamMetricSummary.duplicateReads() / bamMetricSummary.totalReads(),
                QcThresholdRegistry.getCommonThreshold(DUPLICATE_READS));

        putFeature(featuresMap, DUAL_STRAND_READS, (double) bamMetricSummary.dualStrandReads() / bamMetricSummary.totalReads(),
                QcThresholdRegistry.getCommonThreshold(DUAL_STRAND_READS));
    }

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamMetricCoverage bamMetricCoverage, SampleType sampleType)
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

            QcThreshold qcThreshold = QcThresholdRegistry.getThreshold(sampleType, coverageAboveXFeature);
            QcStatus qcStatus = qcThreshold.getQcStatus(propBasesAboveCoverage);

            putFeature(featuresMap, coverageAboveXFeature, propBasesAboveCoverage, qcStatus);
        }
    }

    private static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamFlagStats bamFlagStats)
    {
        if(bamFlagStats == null)
            return;

        QcStatus qcStatus = QcThresholdRegistry.getThreshold(ThresholdGroup.COMMON, MAPPED_PROPORTION).getQcStatus(bamFlagStats.mappedProportion());
        putFeature(featuresMap, MAPPED_PROPORTION, bamFlagStats.mappedProportion(), qcStatus);
    }

    @VisibleForTesting
    static List<Feature> createFeatures(BamMetricsData bamMetricsData, SampleType sampleType)
    {
        EnumMap<SummaryTableFeature, Feature> featuresMap = new EnumMap<>(SummaryTableFeature.class);

        putFeatures(featuresMap, bamMetricsData.bamMetricSummary(), sampleType);
        putFeatures(featuresMap, bamMetricsData.bamMetricCoverage(), sampleType);
        putFeatures(featuresMap, bamMetricsData.bamFlagStats());

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

        List<Feature> features = createFeatures(bamMetricsData, sampleType);

        return features;
    }

}
