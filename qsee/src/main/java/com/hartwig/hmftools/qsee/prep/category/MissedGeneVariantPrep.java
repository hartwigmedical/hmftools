package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.jetbrains.annotations.NotNull;

public class MissedGeneVariantPrep implements CategoryPrep
{
    private final QseePrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    private static final String FIELD_GENE = "Gene";

    private static final int MIN_VARIANT_ALLELE_DEPTH = 3;
    private static final double REF_SAMPLE_PURITY = 1;  // Assume no contamination from other samples

    private static final int LOD_VAF_BISECTION_ITERATIONS = 20; // VAF precision 1 / 2^20 = ~1e-6
    private static final double LOD_TARGET_MISSED_VARIANT_LIKELIHOOD = 0.01;

    public MissedGeneVariantPrep(QseePrepConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }
    public PrepCategory category() { return PrepCategory.MISSED_GENE_VARIANT; }

    private List<GeneDepth> loadGeneCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = GeneDepthFile.generateGeneCoverageFilename(baseDir, sampleId);
        return GeneDepthFile.read(filePath);
    }

    private PurplePurity loadPurplePurity(String sampleId) throws IOException
    {
        String baseDir = mConfig.getPurpleDir(sampleId);
        String purityFile = PurplePurity.generateFilename(baseDir, sampleId);
        return PurplePurity.read(purityFile);
    }

    @VisibleForTesting
    static List<String> getReportableGenes(List<DriverGene> driverGenes, SampleType sampleType)
    {
        List<String> reportableGenes = new ArrayList<>();

        for(DriverGene driverGene : driverGenes)
        {
            boolean isRelevantReportableGene = sampleType == SampleType.TUMOR
                    ? driverGene.reportSomatic()
                    : driverGene.reportGermline();

            if(isRelevantReportableGene)
            {
                reportableGenes.add(driverGene.gene());
            }
        }

        return reportableGenes;
    }

    @VisibleForTesting
    static double calcMissedVariantLikelihood(Map<DepthRange,Integer> depthCounts, double expectedVAF)
    {
        int numBasesInRegion = depthCounts.values().stream().mapToInt(Integer::intValue).sum();
        double missedVariantLikelihood = 0;

        for(DepthRange depthRange : depthCounts.keySet())
        {
            int numBasesWithDepth = depthCounts.get(depthRange);

            if(numBasesWithDepth == 0)
                continue;

            int numReads = depthRange.midpoint();

            double probMutationNotObserved;
            if(numReads > 0 && expectedVAF > 0)
            {
                BinomialDistribution distribution = new BinomialDistribution(numReads, expectedVAF);
                probMutationNotObserved = distribution.cumulativeProbability(MIN_VARIANT_ALLELE_DEPTH - 1);
            }
            else
            {
                probMutationNotObserved = 1;
            }

            double propBasesWithDepth = (double) numBasesWithDepth / numBasesInRegion;
            missedVariantLikelihood += propBasesWithDepth * probMutationNotObserved;
        }

        // 'Likelihood' is a misnomer (should be 'probability'), but will leave as is because renaming is non-trivial
        return missedVariantLikelihood;
    }

    @VisibleForTesting
    static double calcPurityAdjustedLimitOfDetectionVAF(Map<DepthRange,Integer> depthCounts, double purity)
    {
        // As VAF (x) increases, missed variant likelihood (y) monotonically decreases (i.e. downward slope)
        // We can therefore use a binary search approach (i.e. bisection) to quickly find 'x' given 'y'

        double minAchievableMissedVariantLikelihood = calcMissedVariantLikelihood(depthCounts, 1);
        if(minAchievableMissedVariantLikelihood > LOD_TARGET_MISSED_VARIANT_LIKELIHOOD)
            return Double.POSITIVE_INFINITY;

        double lowerVAF = 0;
        double upperVAF = 1;

        for(int iteration = 0; iteration < LOD_VAF_BISECTION_ITERATIONS; ++iteration)
        {
            double midVAF = (lowerVAF + upperVAF) / 2;

            if(calcMissedVariantLikelihood(depthCounts, midVAF) <= LOD_TARGET_MISSED_VARIANT_LIKELIHOOD)
                upperVAF = midVAF;
            else
                lowerVAF = midVAF;
        }

        double limitOfDetectionVAF = upperVAF;

        // Lower purity means higher theoretical VAF required to detect variant in tumor
        double adjustedLimitOfDetectionVAF = limitOfDetectionVAF / purity;

        if(adjustedLimitOfDetectionVAF > 1)
            return Double.POSITIVE_INFINITY;

        return adjustedLimitOfDetectionVAF;
    }

    @VisibleForTesting
    static List<Feature> getMissedVariantLikelihoods(List<GeneDepth> geneDepths, List<String> reportableGenes, double purity, Gender gender)
    {
        List<GeneDepth> selectedGeneDepths = geneDepths.stream().filter(x -> reportableGenes.contains(x.Gene)).toList();

        List<Feature> features = new ArrayList<>();
        for(GeneDepth geneDepth : selectedGeneDepths)
        {
            String featureName = MultiFieldStringBuilder.formSingleField(FIELD_GENE, geneDepth.Gene);
            FeatureKey key = new FeatureKey(featureName, FeatureType.MISSED_VARIANT_LIKELIHOOD, SOURCE_TOOL);

            HumanChromosome geneChromosome = HumanChromosome.fromString(geneDepth.Chromosome);
            int alleleCount = geneChromosome.isDiploid(gender) ? 2 : 1;

            // Assume there is 1 mutation on 1 allele at a particular position.
            // With 0.5 purity and 2 alleles (diploid), the prior probability of observing the mutation on a read (i.e. expected VAF) is:
            //   0.5 purity * (1 mutation / 2 alleles) = 0.25
            double expectedVAF = purity * (1.0 / alleleCount);

            double missedVariantLikelihood = calcMissedVariantLikelihood(geneDepth.DepthCounts, expectedVAF);

            Feature feature = new Feature(key, missedVariantLikelihood);
            features.add(feature);
        }

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        List<GeneDepth> geneCoverage = loadGeneCoverage(sampleId, sampleType);
        List<String> reportableGenes = getReportableGenes(mConfig.DriverGenes, sampleType);

        double purity;
        Gender gender;

        if(sampleType == SampleType.TUMOR)
        {
            PurplePurity purplePurity = loadPurplePurity(sampleId);
            purity = purplePurity.Purity;
            gender = purplePurity.Sex;
        }
        else
        {
            purity = REF_SAMPLE_PURITY;
            // Reference samples don't have Purple output - use that of the paired tumor sample to get gender
            String tumorId = mConfig.TumorIds.get(mConfig.ReferenceIds.indexOf(sampleId));
            gender = loadPurplePurity(tumorId).Sex;
        }

        List<Feature> features = getMissedVariantLikelihoods(geneCoverage, reportableGenes, purity, gender);
        return features;
    }
}
