package prep.category;

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

import common.SampleType;
import feature.Feature;
import feature.FeatureType;
import feature.SourceTool;
import prep.CategoryPrep;
import prep.PrepConfig;

public class SummaryTablePrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    public SummaryTablePrep(PrepConfig config)
    {
        mConfig = config;
    }

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

    private static List<Feature> getGeneralStats(PurityContext purplePurityContext, AmberQC amberQC)
    {
        return List.of(
                new Feature(FeatureType.GENERAL, "Purity", purplePurityContext.bestFit().purity(), SourceTool.PURPLE),
                new Feature(FeatureType.GENERAL, "Ploidy", purplePurityContext.bestFit().ploidy(), SourceTool.PURPLE),
                new Feature(FeatureType.GENERAL, "TINC level", purplePurityContext.qc().tincLevel(), SourceTool.PURPLE),
                new Feature(FeatureType.GENERAL, "Deleted genes", purplePurityContext.qc().deletedGenes(), SourceTool.PURPLE),
                new Feature(FeatureType.GENERAL, "Unsupported CN segments", purplePurityContext.qc().unsupportedCopyNumberSegments(), SourceTool.PURPLE),
                new Feature(FeatureType.GENERAL, "LOH percent", purplePurityContext.qc().lohPercent(), SourceTool.PURPLE),
                new Feature(FeatureType.GENERAL, "Chimerism percent", purplePurityContext.qc().chimerismPercentage(), SourceTool.PURPLE),
                new Feature(FeatureType.GENERAL, "Contamination", purplePurityContext.qc().contamination(), SourceTool.PURPLE),
                new Feature(FeatureType.GENERAL, "Consanguinity prop.", amberQC.consanguinityProportion(), SourceTool.PURPLE)
        );
    }

    private static List<Feature> getTmbStats(PurityContext purplePurityContext)
    {
        return List.of(
                new Feature(FeatureType.MUTATIONAL_BURDEN, "Small variants per MB", purplePurityContext.tumorMutationalBurdenPerMb(), SourceTool.PURPLE),
                new Feature(FeatureType.MUTATIONAL_BURDEN, "MS indels per per MB", purplePurityContext.microsatelliteIndelsPerMb(), SourceTool.PURPLE),
                new Feature(FeatureType.MUTATIONAL_BURDEN, "Structural variants", purplePurityContext.svTumorMutationalBurden(), SourceTool.PURPLE)
        );

    }

    private static final int[] MIN_COVERAGE_THRESHOLDS = { 10, 30, 100, 250 };

    private static List<Feature> getCoverageStats(BamMetricSummary bamMetricSummary, BamMetricCoverage bamMetricCoverage)
    {
        List<Feature> features = new ArrayList<>();

        features.add(new Feature(FeatureType.COVERAGE_STATS, "Mean coverage", bamMetricSummary.meanCoverage(), SourceTool.PURPLE));

        for(int coverageThreshold : MIN_COVERAGE_THRESHOLDS)
        {
            features.add(calcPropBasesWithMinCoverage(bamMetricCoverage, coverageThreshold));
        }

        features.add(new Feature(FeatureType.COVERAGE_STATS, "Low map qual percent", bamMetricSummary.lowMapQualPercent(), SourceTool.BAM_METRICS));
        features.add(new Feature(FeatureType.COVERAGE_STATS, "Low base qual percent", bamMetricSummary.lowBaseQualPercent(), SourceTool.BAM_METRICS));

        return features;
    }


    private static Feature calcPropBasesWithMinCoverage(BamMetricCoverage bamMetricCoverage, int coverageThreshold)
    {
        List<ValueFrequency> coverageBaseCounts = bamMetricCoverage.Coverage;

        long totalBases = coverageBaseCounts.stream().mapToLong(x -> x.Count).sum();

        long basesAboveCoverageThres = coverageBaseCounts.stream()
                .filter(x -> x.Value >= coverageThreshold)
                .mapToLong(x -> x.Count)
                .sum();

        double propAboveCoverageThres = (double) basesAboveCoverageThres / totalBases;

        return new Feature(
                FeatureType.COVERAGE_STATS,
                String.format("Coverage ≥ %s", coverageThreshold),
                propAboveCoverageThres,
                SourceTool.BAM_METRICS
        );
    }

    private static List<Feature> getReadStats(BamMetricSummary bamMetricSummary)
    {
        return List.of(
                new Feature(FeatureType.READ_STATS, "Duplicate reads rate", (double) bamMetricSummary.duplicateReads() / bamMetricSummary.totalReads(), SourceTool.BAM_METRICS),
                new Feature(FeatureType.READ_STATS, "Dual-strand reads rate", (double) bamMetricSummary.dualStrandReads() / bamMetricSummary.totalReads(), SourceTool.BAM_METRICS)
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
