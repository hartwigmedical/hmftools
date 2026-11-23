package com.hartwig.hmftools.qsee.feature.category;

import static com.hartwig.hmftools.qsee.coverage.LimitOfDetectionCalcs.MIN_ALT_READS;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.qsee.QseeConfig;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.common.ReportableGenes;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.coverage.ExonExclusions;
import com.hartwig.hmftools.qsee.coverage.LimitOfDetectionCalcs;
import com.hartwig.hmftools.qsee.feature.CategoryPrep;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

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

    private List<ExonDepth> loadExonCoverage(String sampleId, SampleType sampleType) throws IOException
    {
        String baseDir = mConfig.getBamMetricsDir(sampleId, sampleType);
        String filePath = GeneDepthFile.generateExonCoverageFilename(baseDir, sampleId);
        return GeneDepthFile.readExonCoverage(filePath);
    }

    private PurplePurity loadPurplePurity(String sampleId) throws IOException
    {
        String baseDir = mConfig.getPurpleDir(sampleId);
        String purityFile = PurplePurity.generateFilename(baseDir, sampleId);
        return PurplePurity.read(purityFile);
    }

    @VisibleForTesting
    static List<Feature> getMissedVariantLikelihoods(
            List<ExonDepth> exonDepths,
            List<String> reportableGenes,
            ExonExclusions exonExclusions,
            double purity, Gender gender
    ){
        // Gene depth counts are the sum of the non-excluded exon depth counts
        Map<String,Map<DepthRange,Integer>> geneDepthCounts = new LinkedHashMap<>();
        Map<String,String> geneChromosomes = new HashMap<>();

        for(ExonDepth exonDepth : exonDepths)
        {
            if(!reportableGenes.contains(exonDepth.Gene))
                continue;

            if(exonExclusions.isExcluded(exonDepth.Gene, exonDepth.ExonRank))
                continue;

            Map<DepthRange,Integer> depthCounts = geneDepthCounts.computeIfAbsent(exonDepth.Gene, k -> new LinkedHashMap<>());
            exonDepth.DepthCounts.forEach((depthRange, count) -> depthCounts.merge(depthRange, count, Integer::sum));

            geneChromosomes.put(exonDepth.Gene, exonDepth.Chromosome);
        }

        List<Feature> features = new ArrayList<>();
        for(Map.Entry<String,Map<DepthRange,Integer>> entry : geneDepthCounts.entrySet())
        {
            String gene = entry.getKey();

            String featureName = MultiFieldStringBuilder.formSingleField(FIELD_GENE, gene);
            FeatureKey key = new FeatureKey(featureName, FeatureType.MISSED_VARIANT_LIKELIHOOD, SOURCE_TOOL);

            HumanChromosome geneChromosome = HumanChromosome.fromString(geneChromosomes.get(gene));
            double theoreticalVaf = LimitOfDetectionCalcs.calcTheoreticalVafMutationOnOneAllele(purity, geneChromosome, gender);

            double missedVariantLikelihood = LimitOfDetectionCalcs.calcMissedVariantLikelihood(entry.getValue(), theoreticalVaf, MIN_ALT_READS);

            Feature feature = new Feature(key, missedVariantLikelihood);
            features.add(feature);
        }

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, @NotNull SampleType sampleType) throws IOException
    {
        List<ExonDepth> exonCoverage = loadExonCoverage(sampleId, sampleType);
        List<String> reportableGenes = ReportableGenes.forSampleType(mConfig.DriverGenes, sampleType);

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

        List<Feature> features = getMissedVariantLikelihoods(exonCoverage, reportableGenes, mConfig.ExcludedExons, purity, gender);
        return features;
    }
}
