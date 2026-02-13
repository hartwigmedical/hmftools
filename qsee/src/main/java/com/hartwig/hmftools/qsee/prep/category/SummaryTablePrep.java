package com.hartwig.hmftools.qsee.prep.category;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeConstants.MB_PER_GENOME;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.PurpleQCFile;

import org.jetbrains.annotations.NotNull;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CategoryPrepTask;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;

public class SummaryTablePrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.MULTIPLE;

    public SummaryTablePrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private PurityContext loadPurplePurity(String sampleId, List<String> missingInputPaths)
    {
        String baseDir = mConfig.getPurpleDir(sampleId);
        String purityFile = PurplePurity.generateFilename(baseDir, sampleId);
        String qcFile = PurpleQCFile.generateFilename(baseDir, sampleId);

        try
        {
            return PurityContextFile.readWithQC(qcFile, purityFile);
        }
        catch(IOException e)
        {
            missingInputPaths.add(purityFile);
            missingInputPaths.add(qcFile);
            return null;
        }
    }

    private BamMetricSummary loadBamMetricSummary(String sampleId, SampleType sampleType, List<String> missingInputPaths)
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricSummary.generateFilename(baseDir, sampleId);

        try
        {
            return BamMetricSummary.read(filePath);
        }
        catch(IOException e)
        {
            missingInputPaths.add(filePath);
            return null;
        }
    }

    private BamMetricCoverage loadBamMetricCoverage(String sampleId, SampleType sampleType, List<String> missingInputPaths)
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricCoverage.generateFilename(baseDir, sampleId);

        try
        {
            return BamMetricCoverage.read(filePath);
        }
        catch(IOException e)
        {
            missingInputPaths.add(filePath);
            return null;
        }
    }

    @VisibleForTesting
    private static void putFeature(EnumMap<SummaryTableFeature, Feature> featuresMap, SummaryTableFeature summaryTableFeature, double value)
    {
        Feature feature = new Feature(summaryTableFeature.key(), value);
        featuresMap.put(summaryTableFeature, feature);
    }

    @VisibleForTesting
    static void putFeatures(PurityContext purityContext, EnumMap<SummaryTableFeature, Feature> featuresMap)
    {
        if(purityContext == null)
            return;

        putFeature(featuresMap, SummaryTableFeature.PURITY, purityContext.bestFit().purity());
        putFeature(featuresMap, SummaryTableFeature.PLOIDY, purityContext.bestFit().ploidy());
        putFeature(featuresMap, SummaryTableFeature.TINC, purityContext.qc().tincLevel());
        putFeature(featuresMap, SummaryTableFeature.DELETED_GENES, purityContext.qc().deletedGenes());
        putFeature(featuresMap, SummaryTableFeature.UNSUPPORTED_CN_SEGMENTS, purityContext.qc().unsupportedCopyNumberSegments());
        putFeature(featuresMap, SummaryTableFeature.LOH_PERCENT, purityContext.qc().lohPercent());
        putFeature(featuresMap, SummaryTableFeature.CONTAMINATION, purityContext.qc().contamination());
        putFeature(featuresMap, SummaryTableFeature.TMB_SMALL_VARIANTS, purityContext.tumorMutationalBurdenPerMb());
        putFeature(featuresMap, SummaryTableFeature.TMB_MS_INDELS, purityContext.microsatelliteIndelsPerMb());
        putFeature(featuresMap, SummaryTableFeature.TMB_STRUCTURAL_VARIANTS, purityContext.svTumorMutationalBurden() / MB_PER_GENOME);
    }

    @VisibleForTesting
    static void putFeatures(BamMetricSummary bamMetricSummary, EnumMap<SummaryTableFeature, Feature> featuresMap)
    {
        if(bamMetricSummary == null)
            return;

        putFeature(featuresMap, SummaryTableFeature.MEAN_COVERAGE, bamMetricSummary.meanCoverage());
        putFeature(featuresMap, SummaryTableFeature.LOW_MAP_QUAL, bamMetricSummary.lowMapQualPercent());
        putFeature(featuresMap, SummaryTableFeature.LOW_BASE_QUAL, bamMetricSummary.lowBaseQualPercent());
        putFeature(featuresMap, SummaryTableFeature.DUPLICATE_READS, (double) bamMetricSummary.duplicateReads() / bamMetricSummary.totalReads());
        putFeature(featuresMap, SummaryTableFeature.DUAL_STRAND_READS, (double) bamMetricSummary.dualStrandReads() / bamMetricSummary.totalReads());
    }

    @VisibleForTesting
    static void putFeatures(BamMetricCoverage bamMetricCoverage, EnumMap<SummaryTableFeature, Feature> featuresMap)
    {
        if(bamMetricCoverage == null)
            return;

        putFeature(featuresMap, SummaryTableFeature.MIN_COVERAGE_10, calcPropBasesWithMinCoverage(bamMetricCoverage, 10));
        putFeature(featuresMap, SummaryTableFeature.MIN_COVERAGE_20, calcPropBasesWithMinCoverage(bamMetricCoverage, 20));
        putFeature(featuresMap, SummaryTableFeature.MIN_COVERAGE_30, calcPropBasesWithMinCoverage(bamMetricCoverage, 30));
        putFeature(featuresMap, SummaryTableFeature.MIN_COVERAGE_60, calcPropBasesWithMinCoverage(bamMetricCoverage, 60));
        putFeature(featuresMap, SummaryTableFeature.MIN_COVERAGE_100, calcPropBasesWithMinCoverage(bamMetricCoverage, 100));
        putFeature(featuresMap, SummaryTableFeature.MIN_COVERAGE_250, calcPropBasesWithMinCoverage(bamMetricCoverage, 250));
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

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        EnumMap<SummaryTableFeature, Feature> featuresMap = new EnumMap<>(SummaryTableFeature.class);
        List<String> missingInputPaths = new ArrayList<>();

        if(sampleType == SampleType.TUMOR)
        {
            PurityContext purityContext = loadPurplePurity(sampleId, missingInputPaths);
            putFeatures(purityContext, featuresMap);
        }

        BamMetricSummary bamMetricSummary = loadBamMetricSummary(sampleId, sampleType, missingInputPaths);
        putFeatures(bamMetricSummary, featuresMap);

        BamMetricCoverage bamMetricCoverage = loadBamMetricCoverage(sampleId, sampleType, missingInputPaths);
        putFeatures(bamMetricCoverage, featuresMap);

        if(!missingInputPaths.isEmpty())
        {
            StringJoiner toolsMissingInput = new StringJoiner(", ");
            for(String path : missingInputPaths)
            {
                String basename = new File(path).getName();
                toolsMissingInput.add(basename);
            }

            CategoryPrepTask.missingInputFilesError(
                    mConfig.AllowMissingInput, this, sampleType, sampleId, toolsMissingInput.toString());
        }

        return featuresMap.values().stream().toList();
    }
}
