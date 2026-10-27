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

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;

public class MissedGeneVariantPrep implements CategoryPrep
{
    private final QseePrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    private static final String FIELD_GENE = "Gene";

    private static final int MIN_VARIANT_ALLELE_DEPTH = 3;
    private static final double REF_SAMPLE_PURITY = 1;  // Assume no contamination from other samples

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

    private static double calcMissedVariantLikelihood(GeneDepth geneDepth, double purity, Gender gender)
    {
        int alleleCount = HumanChromosome.fromString(geneDepth.Chromosome).isDiploid(gender) ? 2 : 1;

        Map<DepthRange,Integer> depthCounts = geneDepth.DepthCounts;
        int totalBaseCount = depthCounts.values().stream().mapToInt(Integer::intValue).sum();
        double totalLikelihood = 0;

        for(DepthRange depthRange : depthCounts.keySet())
        {
            int baseCount = depthCounts.get(depthRange);

            if(baseCount == 0)
                continue;

            double meanAlleleDepth = (purity * depthRange.midpoint()) / alleleCount;

            double likelihoodOfMissing;
            if(meanAlleleDepth > 0)
            {
                PoissonDistribution distribution = new PoissonDistribution(meanAlleleDepth);
                likelihoodOfMissing = distribution.cumulativeProbability(MIN_VARIANT_ALLELE_DEPTH - 1);
            }
            else
            {
                likelihoodOfMissing = 1;
            }

            double baseProportion = (double) baseCount / totalBaseCount;
            totalLikelihood += baseProportion * likelihoodOfMissing;
        }

        return totalLikelihood;
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

            double missedVariantLikelihood = calcMissedVariantLikelihood(geneDepth, purity, gender);

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
