package com.hartwig.hmftools.qsee.prep.category;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeConstants.MB_PER_GENOME;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.*;

import java.io.File;
import java.io.IOException;
import java.util.EnumMap;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.common.purple.PurityContext;

import org.jetbrains.annotations.NotNull;

import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.PlotMetadata;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CategoryPrepTask;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableInputData;
import com.hartwig.hmftools.qsee.status.ComparisonOperator;
import com.hartwig.hmftools.qsee.status.QcStatusType;
import com.hartwig.hmftools.qsee.status.QcThresholdRegistry;
import com.hartwig.hmftools.qsee.status.QcStatus;
import com.hartwig.hmftools.qsee.status.QcThreshold;
import com.hartwig.hmftools.qsee.status.ThresholdGroup;

public class SummaryTablePrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.MULTIPLE;

    public SummaryTablePrep(CommonPrepConfig config) { mConfig = config; }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        SummaryTableInputData inputData = new SummaryTableInputData(mConfig, sampleId, sampleType);

        EnumMap<SummaryTableFeature, Feature> featuresMap = new EnumMap<>(SummaryTableFeature.class);

        putFeatures(featuresMap, inputData.purityContext());
        putFeatures(featuresMap, inputData.bamMetricSummary(), sampleType);
        putFeatures(featuresMap, inputData.bamMetricCoverage(), sampleType);
        putFeatures(featuresMap, inputData.bamFlagStats());

        printMissingInputFiles(inputData, sampleType, sampleId);

        return featuresMap.values().stream().toList();
    }

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, PurityContext purityContext)
    {
        if(purityContext == null)
            return;

        EnumMap<PurpleQCStatus, QcStatus> qcStatuses = qcStatusFrom(purityContext);

        putFeature(featuresMap, PURITY, purityContext.bestFit().purity(),
                qcStatuses.get(PurpleQCStatus.WARN_LOW_PURITY));

        putFeature(featuresMap, PLOIDY, purityContext.bestFit().ploidy());

        putFeature(featuresMap, TINC, purityContext.qc().tincLevel(),
                getTincQcStatus(qcStatuses));

        putFeature(featuresMap, DELETED_GENES, purityContext.qc().deletedGenes(),
                qcStatuses.get(PurpleQCStatus.WARN_DELETED_GENES));

        putFeature(featuresMap, UNSUPPORTED_CN_SEGMENTS, purityContext.qc().unsupportedCopyNumberSegments(),
                qcStatuses.get(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE));

        putFeature(featuresMap, LOH_PERCENT, purityContext.qc().lohPercent());

        putFeature(featuresMap, CONTAMINATION, purityContext.qc().contamination(),
                qcStatuses.get(PurpleQCStatus.FAIL_CONTAMINATION));

        putFeature(featuresMap, TMB_SMALL_VARIANTS, purityContext.tumorMutationalBurdenPerMb());
        putFeature(featuresMap, TMB_MS_INDELS, purityContext.microsatelliteIndelsPerMb());
        putFeature(featuresMap, TMB_STRUCTURAL_VARIANTS, purityContext.svTumorMutationalBurden() / MB_PER_GENOME);
    }

    @VisibleForTesting
    static EnumMap<PurpleQCStatus, QcStatus> qcStatusFrom(PurityContext purityContext)
    {
        EnumMap<PurpleQCStatus, QcStatus> qcStatuses = new EnumMap<>(PurpleQCStatus.class);

        for(PurpleQCStatus purpleQcStatus : PurpleQCStatus.values())
        {
            boolean purpleQcStatusExists = purityContext.qc().status().contains(purpleQcStatus);

            QcStatus qcStatus = purpleQcStatusExists
                    ? convertQcStatus(purpleQcStatus)
                    : QcStatus.createEmpty();

            qcStatuses.put(purpleQcStatus, qcStatus);
        }

        return qcStatuses;
    }

    private static QcStatus convertQcStatus(PurpleQCStatus purpleQCStatus)
    {
        return switch(purpleQCStatus)
        {
            case PASS -> QcStatus.createEmpty();

            case WARN_DELETED_GENES ->
                    new QcStatus(QcStatusType.WARN, ComparisonOperator.GREATER_THAN, PurpleQCStatus.MAX_DELETED_GENES);

            case WARN_HIGH_COPY_NUMBER_NOISE ->
                    new QcStatus(QcStatusType.WARN, ComparisonOperator.GREATER_THAN, PurpleQCStatus.MAX_UNSUPPORTED_SEGMENTS);

            case WARN_LOW_PURITY ->
                    new QcStatus(QcStatusType.WARN, ComparisonOperator.LESS_THAN, PurpleQCStatus.MIN_PURITY);

            case WARN_TINC ->
                    new QcStatus(QcStatusType.WARN, ComparisonOperator.GREATER_THAN, PurpleQCStatus.TINC_WARN_LEVEL);

            case FAIL_TINC ->
                    new QcStatus(QcStatusType.FAIL, ComparisonOperator.GREATER_THAN_OR_EQUAL, PurpleQCStatus.TINC_FAIL_LEVEL);

            case FAIL_CONTAMINATION ->
                    new QcStatus(QcStatusType.FAIL, ComparisonOperator.GREATER_THAN, PurpleQCStatus.MAX_CONTAMINATION);

            case FAIL_NO_TUMOR -> new QcStatus(QcStatusType.FAIL, null, Double.NaN);

            case WARN_GENDER_MISMATCH -> new QcStatus(QcStatusType.WARN, null, Double.NaN);

            default ->
                    throw new IllegalArgumentException("QC threshold not defined for PurpleQCStatus(" + purpleQCStatus + ")");
        };
    }

    private static QcStatus getTincQcStatus(EnumMap<PurpleQCStatus, QcStatus> qcStatuses)
    {
        QcStatus tincStatus = qcStatuses.get(PurpleQCStatus.FAIL_NO_TUMOR);

        if(tincStatus == null)
            tincStatus = qcStatuses.get(PurpleQCStatus.WARN_TINC);

        if(tincStatus == null)
            tincStatus = QcStatus.createEmpty();

        return tincStatus;
    }

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamMetricSummary bamMetricSummary, SampleType sampleType)
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

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamFlagStats bamFlagStats)
    {
        if(bamFlagStats == null)
            return;

        QcStatus qcStatus = QcThresholdRegistry.getThreshold(ThresholdGroup.COMMON, MAPPED_PROPORTION).getQcStatus(bamFlagStats.mappedProportion());
        putFeature(featuresMap, MAPPED_PROPORTION, bamFlagStats.mappedProportion(), qcStatus);
    }

    private static void putFeature(
            EnumMap<SummaryTableFeature, Feature> featuresMap,
            SummaryTableFeature summaryTableFeature,
            double value,
            QcStatus qcStatus
    ){
        FeatureKey key = new FeatureKey(summaryTableFeature.toString(), FeatureType.SUMMARY_TABLE, SOURCE_TOOL);

        PlotMetadata metadata = PlotMetadata.builder()
                .featureGroup(summaryTableFeature.group().humanReadableName())
                .plotLabel(summaryTableFeature.plotLabel())
                .numberFormat(summaryTableFeature.numberFormat())
                .qcStatus(qcStatus)
                .build();

        Feature feature = new Feature(key, value, metadata);
        featuresMap.put(summaryTableFeature, feature);
    }

    private static void putFeature(
            EnumMap<SummaryTableFeature, Feature> featuresMap,
            SummaryTableFeature summaryTableFeature,
            double value,
            QcThreshold qcThreshold
    ){
        putFeature(featuresMap, summaryTableFeature, value, qcThreshold.getQcStatus(value));
    }

    private static void putFeature(
            EnumMap<SummaryTableFeature, Feature> featuresMap,
            SummaryTableFeature summaryTableFeature,
            double value
    ){
        putFeature(featuresMap, summaryTableFeature, value, QcStatus.createEmpty());
    }

    private void printMissingInputFiles(SummaryTableInputData inputData, @NotNull SampleType sampleType, String sampleId)
    {
        if(inputData.missingInputPaths().isEmpty())
            return;

        StringJoiner toolsMissingInput = new StringJoiner(", ");
        for(String path : inputData.missingInputPaths())
        {
            String basename = new File(path).getName();
            toolsMissingInput.add(basename);
        }

        CategoryPrepTask.missingInputFilesError(mConfig.AllowMissingInput, this, sampleType, sampleId, toolsMissingInput.toString());
    }
}
