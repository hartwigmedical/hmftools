package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.amber.AmberQC;
import com.hartwig.hmftools.common.amber.AmberQCFile;
import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.ValueFrequency;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleQCFile;

import org.jetbrains.annotations.NotNull;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;
import com.hartwig.hmftools.qsee.table.SummaryTableFeature;

public class SummaryTablePrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final String NAME = "Summary table";
    private static final SourceTool SOURCE_TOOL = SourceTool.MULTIPLE;

    public SummaryTablePrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public String name() { return NAME; }
    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private AmberQC loadAmberQC(String sampleId) throws IOException
    {
        String baseDir = mConfig.getAmberDir(sampleId);
        String filePath = AmberQCFile.generateFilename(baseDir, sampleId);
        return AmberQCFile.read(filePath);
    }

    private PurityContext loadPurplePurityContext(String sampleId) throws IOException
    {
        String baseDir = mConfig.getPurpleDir(sampleId);
        String purityFile = PurityContextFile.generateFilenameForReading(baseDir, sampleId);
        String qcFile = PurpleQCFile.generateFilename(baseDir, sampleId);
        return PurityContextFile.readWithQC(qcFile, purityFile);
    }

    private BamMetricSummary loadBamMetricSummary(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricSummary.generateFilename(baseDir, sampleId);
        return BamMetricSummary.read(filePath);
    }

    private BamMetricCoverage loadBamMetricCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricCoverage.generateFilename(baseDir, sampleId);
        return BamMetricCoverage.read(filePath);
    }

    private static Feature createFeature(SummaryTableFeature feature, double value)
    {
        return new Feature(feature.key(), value, feature.sourceTool());
    }

    private static List<Feature> getGeneralStats(PurityContext purplePurityContext, AmberQC amberQC)
    {
        return List.of(
                createFeature(SummaryTableFeature.PURITY, purplePurityContext.bestFit().purity()),
                createFeature(SummaryTableFeature.PLOIDY, purplePurityContext.bestFit().ploidy()),
                createFeature(SummaryTableFeature.TINC, purplePurityContext.qc().tincLevel()),
                createFeature(SummaryTableFeature.DELETED_GENES, purplePurityContext.qc().deletedGenes()),
                createFeature(SummaryTableFeature.UNSUPPORTED_CN_SEGMENTS, purplePurityContext.qc().unsupportedCopyNumberSegments()),
                createFeature(SummaryTableFeature.LOH_PERCENT, purplePurityContext.qc().lohPercent()),
                createFeature(SummaryTableFeature.CHIMERISM_PERCENT, purplePurityContext.qc().chimerismPercentage()),
                createFeature(SummaryTableFeature.CONTAMINATION, purplePurityContext.qc().contamination()),
                createFeature(SummaryTableFeature.CONSANGUINITY, amberQC.consanguinityProportion())
        );
    }

    private static List<Feature> getTmbStats(PurityContext purplePurityContext)
    {
        return List.of(
                createFeature(SummaryTableFeature.TMB_SMALL_VARIANTS, purplePurityContext.tumorMutationalBurdenPerMb()),
                createFeature(SummaryTableFeature.TMB_MS_INDELS, purplePurityContext.microsatelliteIndelsPerMb()),
                createFeature(SummaryTableFeature.TMB_STRUCTURAL_VARIANTS, purplePurityContext.svTumorMutationalBurden())
        );
    }

    private static List<Feature> getCoverageStats(BamMetricSummary bamMetricSummary, BamMetricCoverage bamMetricCoverage)
    {
        return List.of(
                createFeature(SummaryTableFeature.MEAN_COVERAGE, bamMetricSummary.meanCoverage()),
                createFeature(SummaryTableFeature.MIN_COVERAGE_10, calcPropBasesWithMinCoverage(bamMetricCoverage, 10)),
                createFeature(SummaryTableFeature.MIN_COVERAGE_30, calcPropBasesWithMinCoverage(bamMetricCoverage, 30)),
                createFeature(SummaryTableFeature.MIN_COVERAGE_100, calcPropBasesWithMinCoverage(bamMetricCoverage, 100)),
                createFeature(SummaryTableFeature.MIN_COVERAGE_250, calcPropBasesWithMinCoverage(bamMetricCoverage, 250)),
                createFeature(SummaryTableFeature.LOW_MAP_QUAL, bamMetricSummary.lowMapQualPercent()),
                createFeature(SummaryTableFeature.LOW_BASE_QUAL, bamMetricSummary.lowBaseQualPercent())
        );
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

    private static List<Feature> getReadStats(BamMetricSummary bamMetricSummary)
    {
        return List.of(
                createFeature(SummaryTableFeature.DUPLICATE_READS, (double) bamMetricSummary.duplicateReads() / bamMetricSummary.totalReads()),
                createFeature(SummaryTableFeature.DUAL_STRAND_READS, (double) bamMetricSummary.dualStrandReads() / bamMetricSummary.totalReads())
        );
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        List<Feature> features = new ArrayList<>();

        if(sampleType == SampleType.TUMOR)
        {
            AmberQC amberQC = loadAmberQC(sampleId);
            PurityContext purplePurityContext = loadPurplePurityContext(sampleId);

            features.addAll(getGeneralStats(purplePurityContext, amberQC));
            features.addAll(getTmbStats(purplePurityContext));
        }

        BamMetricSummary bamMetricSummary = loadBamMetricSummary(sampleId, sampleType);
        BamMetricCoverage bamMetricCoverage = loadBamMetricCoverage(sampleId, sampleType);

        features.addAll(getCoverageStats(bamMetricSummary, bamMetricCoverage));
        features.addAll(getReadStats(bamMetricSummary));

        return features;
    }
}
