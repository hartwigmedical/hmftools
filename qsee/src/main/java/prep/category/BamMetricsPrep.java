package prep.category;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.BamMetricFragmentLength;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.common.metrics.ValueFrequency;

import org.jetbrains.annotations.NotNull;

import common.SampleType;
import feature.Feature;
import feature.FeatureType;
import prep.CategoryPrep;
import prep.PrepConfig;

public class BamMetricsPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    private BamMetricSummary mBamMetricSummary;
    private BamMetricCoverage mBamMetricCoverage;
    private BamMetricFragmentLength mBamMetricFragmentLength;
    private List<GeneDepth> mBamMetricGeneCoverage;

    public BamMetricsPrep(PrepConfig config)
    {
        mConfig = config;
    }

    private void loadBamMetricsFiles(String bamMetricsDir, String sampleId) throws IOException
    {
        String filePath;

        filePath = BamMetricSummary.generateFilename(bamMetricsDir, sampleId);
        mBamMetricSummary = BamMetricSummary.read(filePath);

        filePath = BamMetricCoverage.generateFilename(bamMetricsDir, sampleId);
        mBamMetricCoverage = BamMetricCoverage.read(filePath);

        filePath = BamMetricFragmentLength.generateFilename(bamMetricsDir, sampleId);
        mBamMetricFragmentLength = BamMetricFragmentLength.read(filePath);

        filePath = GeneDepthFile.generateGeneCoverageFilename(bamMetricsDir, sampleId);
        mBamMetricGeneCoverage = GeneDepthFile.read(filePath);
    }

    private static List<Feature> getSummaryStats(BamMetricSummary summary)
    {
        List<Feature> features = new ArrayList<>();

        features.add(new Feature(FeatureType.COVERAGE_STATS, "MeanCoverage", summary.meanCoverage()));
        features.add(new Feature(FeatureType.COVERAGE_STATS, "LowMapQualPercent", summary.lowMapQualPercent()));
        features.add(new Feature(FeatureType.COVERAGE_STATS, "LowBaseQualPercent", summary.lowBaseQualPercent()));

        features.add(new Feature(FeatureType.READ_STATS, "DuplicateReadsRate", (double) summary.duplicateReads() / summary.totalReads()));
        features.add(new Feature(FeatureType.READ_STATS, "DualStrandReadsRate", (double) summary.dualStrandReads() / summary.totalReads()));

        return features;
    }

    private static final int[] MIN_COVERAGE_THRESHOLDS = { 10, 30, 100, 250 };

    private static List<Feature> calcPropBasesWithMinCoverage(List<ValueFrequency> coverageBaseCounts)
    {
        List<Feature> features = new ArrayList<>();

        long totalBases = coverageBaseCounts.stream().mapToLong(x -> x.Count).sum();
        for(int coverageThreshold : MIN_COVERAGE_THRESHOLDS)
        {
            long basesAboveCoverageThres = coverageBaseCounts.stream()
                    .filter(x -> x.Value >= coverageThreshold)
                    .mapToLong(x -> x.Count)
                    .sum();

            double propAboveCoverageThres = (double) basesAboveCoverageThres / totalBases;

            Feature feature = new Feature(
                    FeatureType.COVERAGE_STATS,
                    String.format("Coverage ≥ %s", coverageThreshold),
                    propAboveCoverageThres
            );

            features.add(feature);
        }

        return features;
    }

    private static List<Feature> calcPropBasesWithCoverage(List<ValueFrequency> coverageBaseCounts)
    {
        long totalCount = coverageBaseCounts.stream().mapToLong(x -> x.Count).sum();

        return coverageBaseCounts.stream().map(x -> {
            double propBases = (double) x.Count / totalCount;
            return new Feature(FeatureType.COVERAGE_DISTRIBUTION, String.valueOf(x.Value), propBases);
        }).toList();
    }

    private static List<Feature> calcPropFragmentsWithLength(List<ValueFrequency> fragmentLengthCounts)
    {
        long totalFragments = fragmentLengthCounts.stream().mapToLong(x -> x.Count).sum();

        return fragmentLengthCounts.stream().map(x -> {
            double propBases = (double) x.Count / totalFragments;
            return new Feature(FeatureType.FRAG_LENGTH_DISTRIBUTION, String.valueOf(x.Value), propBases);
        }).toList();
    }

    private static List<Feature> getMissedVariantLikelihoods(List<GeneDepth> geneDepths, List<DriverGene> driverGenes)
    {
        List<String> selectedGenes = driverGenes.stream().filter(DriverGene::reportSomatic).map(DriverGene::gene).toList();
        List<GeneDepth> selectedGeneDepths = geneDepths.stream().filter(x -> selectedGenes.contains(x.Gene)).toList();

        return selectedGeneDepths.stream()
                .map(x -> new Feature(FeatureType.MISSED_VARIANT_LIKELIHOOD, x.Gene, x.MissedVariantLikelihood))
                .toList();
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        String bamMetricsDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        loadBamMetricsFiles(bamMetricsDir, sampleId);

        List<Feature> summaryStats = getSummaryStats(mBamMetricSummary);
        List<Feature> propBasesWithCoverage = calcPropBasesWithCoverage(mBamMetricCoverage.Coverage);
        List<Feature> propBasesWithMinCoverage  = calcPropBasesWithMinCoverage(mBamMetricCoverage.Coverage);
        List<Feature> fragmentLengthDistribution = calcPropFragmentsWithLength(mBamMetricFragmentLength.FragmentLengths);
        List<Feature> missedVariantLikelihoods = getMissedVariantLikelihoods(mBamMetricGeneCoverage, mConfig.DriverGenes);

        List<Feature> features = new ArrayList<>();
        features.addAll(summaryStats);
        features.addAll(propBasesWithCoverage);
        features.addAll(propBasesWithMinCoverage);
        features.addAll(fragmentLengthDistribution);
        features.addAll(missedVariantLikelihoods);

        return features;
    }
}
