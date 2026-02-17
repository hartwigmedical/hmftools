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

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CategoryPrepTask;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableInputData;

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
        putFeatures(featuresMap, inputData.bamMetricSummary());
        putFeatures(featuresMap, inputData.bamMetricCoverage());
        putFeatures(featuresMap, inputData.bamFlagStats());

        printMissingInputFiles(inputData, sampleType, sampleId);

        return featuresMap.values().stream().toList();
    }

    @VisibleForTesting
    private static void putFeature(EnumMap<SummaryTableFeature, Feature> featuresMap, SummaryTableFeature summaryTableFeature, double value)
    {
        Feature feature = new Feature(summaryTableFeature.key(), value);
        featuresMap.put(summaryTableFeature, feature);
    }

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, PurityContext purityContext)
    {
        if(purityContext == null)
            return;

        putFeature(featuresMap, PURITY, purityContext.bestFit().purity());
        putFeature(featuresMap, PLOIDY, purityContext.bestFit().ploidy());
        putFeature(featuresMap, TINC, purityContext.qc().tincLevel());
        putFeature(featuresMap, DELETED_GENES, purityContext.qc().deletedGenes());
        putFeature(featuresMap, UNSUPPORTED_CN_SEGMENTS, purityContext.qc().unsupportedCopyNumberSegments());
        putFeature(featuresMap, LOH_PERCENT, purityContext.qc().lohPercent());
        putFeature(featuresMap, CONTAMINATION, purityContext.qc().contamination());
        putFeature(featuresMap, TMB_SMALL_VARIANTS, purityContext.tumorMutationalBurdenPerMb());
        putFeature(featuresMap, TMB_MS_INDELS, purityContext.microsatelliteIndelsPerMb());
        putFeature(featuresMap, TMB_STRUCTURAL_VARIANTS, purityContext.svTumorMutationalBurden() / MB_PER_GENOME);
    }

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamMetricSummary bamMetricSummary)
    {
        if(bamMetricSummary == null)
            return;

        putFeature(featuresMap, MEAN_COVERAGE, bamMetricSummary.meanCoverage());
        putFeature(featuresMap, LOW_MAP_QUAL, bamMetricSummary.lowMapQualPercent());
        putFeature(featuresMap, LOW_BASE_QUAL, bamMetricSummary.lowBaseQualPercent());
        putFeature(featuresMap, DUPLICATE_READS, (double) bamMetricSummary.duplicateReads() / bamMetricSummary.totalReads());
        putFeature(featuresMap, DUAL_STRAND_READS, (double) bamMetricSummary.dualStrandReads() / bamMetricSummary.totalReads());
    }

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamMetricCoverage bamMetricCoverage)
    {
        if(bamMetricCoverage == null)
            return;

        putFeature(featuresMap, MIN_COVERAGE_10, calcPropBasesWithMinCoverage(bamMetricCoverage, 10));
        putFeature(featuresMap, MIN_COVERAGE_20, calcPropBasesWithMinCoverage(bamMetricCoverage, 20));
        putFeature(featuresMap, MIN_COVERAGE_30, calcPropBasesWithMinCoverage(bamMetricCoverage, 30));
        putFeature(featuresMap, MIN_COVERAGE_60, calcPropBasesWithMinCoverage(bamMetricCoverage, 60));
        putFeature(featuresMap, MIN_COVERAGE_100, calcPropBasesWithMinCoverage(bamMetricCoverage, 100));
        putFeature(featuresMap, MIN_COVERAGE_250, calcPropBasesWithMinCoverage(bamMetricCoverage, 250));
    }

    private static double calcPropBasesWithMinCoverage(BamMetricCoverage bamMetricCoverage, int coverageThreshold)
    {
        List<ValueFrequency> coverageBaseCounts = bamMetricCoverage.Coverage;

        long totalBases = coverageBaseCounts.stream().mapToLong(x -> x.Count).sum();

        long basesAboveCoverageThres = coverageBaseCounts.stream()
                .filter(x -> x.Value >= coverageThreshold)
                .mapToLong(x -> x.Count)
                .sum();

        return (double) basesAboveCoverageThres / totalBases;
    }

    @VisibleForTesting
    static void putFeatures(EnumMap<SummaryTableFeature, Feature> featuresMap, BamFlagStats bamFlagStats)
    {
        if(bamFlagStats == null)
            return;

        putFeature(featuresMap, MAPPED_PROPORTION, bamFlagStats.mappedProportion());
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
