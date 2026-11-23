package com.hartwig.hmftools.qsee.feature.category;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.qsee.QseeConfig;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.CategoryPrep;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.lod.LimitOfDetectionCalcs;

import org.jetbrains.annotations.NotNull;

public class MissedGeneVariantPrep implements CategoryPrep
{
    private final QseeConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.BAM_METRICS;

    private static final String FIELD_GENE = "Gene";

    public MissedGeneVariantPrep(QseeConfig config)
    {
        mConfig = config;
    }

    public SourceTool sourceTool() { return SOURCE_TOOL; }
    public PrepCategory category() { return PrepCategory.MISSED_GENE_VARIANT; }

    private List<GeneDepth> loadGeneCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = GeneDepthFile.generateGeneCoverageFilename(baseDir, sampleId);
        return GeneDepthFile.readGeneCoverage(filePath);
    }

    private PurplePurity loadPurplePurity(String sampleId) throws IOException
    {
        String baseDir = mConfig.getPurpleDir(sampleId);
        String purityFile = PurplePurity.generateFilename(baseDir, sampleId);
        return PurplePurity.read(purityFile);
    }

    public static List<String> getReportableGenes(List<DriverGene> driverGenes, SampleType sampleType)
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
    static List<Feature> getMissedVariantLikelihoods(List<GeneDepth> geneDepths, List<String> reportableGenes, double purity, Gender gender)
    {
        List<GeneDepth> selectedGeneDepths = geneDepths.stream().filter(x -> reportableGenes.contains(x.Gene)).toList();

        List<Feature> features = new ArrayList<>();
        for(GeneDepth geneDepth : selectedGeneDepths)
        {
            String featureName = MultiFieldStringBuilder.formSingleField(FIELD_GENE, geneDepth.Gene);
            FeatureKey key = new FeatureKey(featureName, FeatureType.MISSED_VARIANT_LIKELIHOOD, SOURCE_TOOL);

            HumanChromosome geneChromosome = HumanChromosome.fromString(geneDepth.Chromosome);
            double theoreticalVaf = LimitOfDetectionCalcs.calcTheoreticalVafMutationOnOneAllele(purity, geneChromosome, gender);

            double missedVariantLikelihood = LimitOfDetectionCalcs.calcMissedVariantLikelihood(geneDepth.DepthCounts, theoreticalVaf);

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
            purity = LimitOfDetectionCalcs.REF_SAMPLE_PURITY;
            // Reference samples don't have Purple output - use that of the paired tumor sample to get gender
            String tumorId = mConfig.TumorIds.get(mConfig.ReferenceIds.indexOf(sampleId));
            gender = loadPurplePurity(tumorId).Sex;
        }

        List<Feature> features = getMissedVariantLikelihoods(geneCoverage, reportableGenes, purity, gender);
        return features;
    }
}
