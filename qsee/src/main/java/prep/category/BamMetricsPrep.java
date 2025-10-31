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

    public BamMetricsPrep(PrepConfig config)
    {
        mConfig = config;
    }

    private BamMetricCoverage loadCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricCoverage.generateFilename(baseDir, sampleId);
        return BamMetricCoverage.read(filePath);
    }

    private BamMetricFragmentLength loadFragmentLengths(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = BamMetricFragmentLength.generateFilename(baseDir, sampleId);
        return BamMetricFragmentLength.read(filePath);
    }

    private List<GeneDepth> loadGeneCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = GeneDepthFile.generateGeneCoverageFilename(baseDir, sampleId);
        return GeneDepthFile.read(filePath);
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
        BamMetricCoverage coverage = loadCoverage(sampleId, sampleType);
        List<Feature> propBasesWithCoverage = calcPropBasesWithCoverage(coverage.Coverage);

        BamMetricFragmentLength fragmentLengths = loadFragmentLengths(sampleId, sampleType);
        List<Feature> fragmentLengthDistribution = calcPropFragmentsWithLength(fragmentLengths.FragmentLengths);

        List<GeneDepth> geneCoverage = loadGeneCoverage(sampleId, sampleType);
        List<Feature> missedVariantLikelihoods = getMissedVariantLikelihoods(geneCoverage, mConfig.DriverGenes);

        List<Feature> features = new ArrayList<>();
        features.addAll(propBasesWithCoverage);
        features.addAll(fragmentLengthDistribution);
        features.addAll(missedVariantLikelihoods);

        return features;
    }
}
