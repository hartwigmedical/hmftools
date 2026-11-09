package com.hartwig.hmftools.qsee.prep.category;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.driver.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

import org.junit.Test;

public class MissedGeneVariantPrepTest
{
    private static final String GERMLINE_GENE = "GERMLINE_GENE";
    private static final String GERMLINE_CHR_X_GENE = "GERMLINE_CHR_X_GENE";
    private static final String SOMATIC_GENE = "SOMATIC_GENE";
    private static final String SOMATIC_CHR_X_GENE = "SOMATIC_CHR_X_GENE";
    private static final String SOMATIC_NO_COVERAGE_GENE = "SOMATIC_NO_COVERAGE_GENE";

    private static final Map<DepthRange,Integer> MIXED_DEPTH_COUNTS = Map.of(
            new DepthRange(0, 0), 0,
            new DepthRange(15, 15), 100,
            new DepthRange(20, 20), 100,
            new DepthRange(40, 49), 800
    );

    private static final Map<DepthRange,Integer> NO_COVERAGE_DEPTH_COUNTS = Map.of(
            new DepthRange(0, 0), 1
    );

    private static final List<DriverGene> DRIVER_GENES = List.of(
            createDriverGene(GERMLINE_GENE, false, true),
            createDriverGene(GERMLINE_CHR_X_GENE, false, true),
            createDriverGene(SOMATIC_GENE, true, false),
            createDriverGene(SOMATIC_CHR_X_GENE, true, false),
            createDriverGene(SOMATIC_NO_COVERAGE_GENE, true, false)
    );

    private static final List<GeneDepth> GENE_DEPTHS = List.of(
            new GeneDepth(GERMLINE_GENE, "1", 100, 200, MIXED_DEPTH_COUNTS),
            new GeneDepth(GERMLINE_CHR_X_GENE, "X", 100, 200, MIXED_DEPTH_COUNTS),
            new GeneDepth(SOMATIC_GENE, "1", 300, 400, MIXED_DEPTH_COUNTS),
            new GeneDepth(SOMATIC_CHR_X_GENE, "X", 300, 400, MIXED_DEPTH_COUNTS),
            new GeneDepth(SOMATIC_NO_COVERAGE_GENE, "1", 500, 600, NO_COVERAGE_DEPTH_COUNTS)
    );

    @Test
    public void canGetGermlineMissedVariantLikelihoods()
    {
        List<String> reportableGenes = MissedGeneVariantPrep.getReportableGenes(DRIVER_GENES, SampleType.NORMAL);
        List<Feature> actualFeatures = MissedGeneVariantPrep.getMissedVariantLikelihoods(GENE_DEPTHS, reportableGenes, 1, Gender.FEMALE);

        List<Feature> expectedFeatures = List.of(
                createExpectedFeature(GERMLINE_GENE, 0.000389385),
                createExpectedFeature(GERMLINE_CHR_X_GENE, 0.000389385) // Diploid chrX in females
        );

        assertFeaturesEqual(expectedFeatures, actualFeatures);
    }

    @Test
    public void canGetSomaticMissedVariantLikelihoods()
    {
        List<String> reportableGenes = MissedGeneVariantPrep.getReportableGenes(DRIVER_GENES, SampleType.TUMOR);
        List<Feature> actualFeatures = MissedGeneVariantPrep.getMissedVariantLikelihoods(GENE_DEPTHS, reportableGenes, 0.5, Gender.MALE);

        List<Feature> expectedFeatures = List.of(
                createExpectedFeature(SOMATIC_GENE, 0.032975381),
                createExpectedFeature(SOMATIC_CHR_X_GENE, 0.000389385), // Haploid chrX in males doubles depth per allele
                createExpectedFeature(SOMATIC_NO_COVERAGE_GENE, 1)
        );

        assertFeaturesEqual(expectedFeatures, actualFeatures);
    }

    private static void assertFeaturesEqual(List<Feature> expectedFeatures, List<Feature> actualFeatures)
    {
        assertEquals(expectedFeatures.size(), actualFeatures.size());

        for(int i = 0; i < expectedFeatures.size(); i++)
        {
            Feature actualFeature = actualFeatures.get(i);
            Feature expectedFeature = expectedFeatures.get(i);
            assertEquals(expectedFeature.key().name(), actualFeature.key().name());
            assertEquals(expectedFeature.value(), actualFeature.value(), 1e-9);
        }
    }

    private static DriverGene createDriverGene(String geneName, boolean reportSomatic, boolean reportGermline)
    {
        DriverGeneGermlineReporting germlineReporting = reportGermline
                ? DriverGeneGermlineReporting.ANY
                : DriverGeneGermlineReporting.NONE;

        return ImmutableDriverGene.builder()
                .gene(geneName)
                .reportMissenseAndInframe(reportSomatic)
                .reportNonsenseAndFrameshift(reportSomatic)
                .reportSplice(reportSomatic)
                .reportDeletion(false)
                .reportHetDeletion(false)
                .reportDisruption(false)
                .hetDeletionThreshold(0.0)
                .reportLoh(false)
                .reportAmplification(false)
                .amplificationRatio(0.0)
                .reportSomaticHotspot(reportSomatic)
                .likelihoodType(DriverCategory.TSG)
                .reportGermlineVariant(germlineReporting)
                .reportGermlineHotspot(germlineReporting)
                .reportGermlineDisruption(germlineReporting)
                .reportGermlineDeletion(germlineReporting)
                .reportGermlineAmplification(false)
                .reportPGX(false)
                .reportHighExpression(false)
                .reportLowExpression(false)
                .reportNovelSpliceJunction(false)
                .build();
    }

    private static Feature createExpectedFeature(String geneName, double missedVariantLikelihood)
    {
        return new Feature("Gene:" + geneName, missedVariantLikelihood, FeatureType.MISSED_VARIANT_LIKELIHOOD, SourceTool.BAM_METRICS, null);
    }
}
