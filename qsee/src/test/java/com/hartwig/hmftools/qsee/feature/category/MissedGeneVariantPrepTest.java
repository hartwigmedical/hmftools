package com.hartwig.hmftools.qsee.feature.category;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.driver.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.metrics.ExonDepth;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.qsee.common.ReportableGenes;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.coverage.ExonExclusions;
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
    private static final String SOMATIC_EXCLUDED_GENE = "SOMATIC_EXCLUDED_GENE";

    private static final Map<DepthRange,Integer> MIXED_DEPTH_COUNTS = Map.of(
            new DepthRange(0, 0), 0,
            new DepthRange(15, 15), 100,
            new DepthRange(20, 20), 100,
            new DepthRange(40, 49), 800
    );

    // Sums to MIXED_DEPTH_COUNTS
    private static final Map<DepthRange,Integer> MIXED_DEPTH_COUNTS_PART_1 = Map.of(
            new DepthRange(0, 0), 0,
            new DepthRange(15, 15), 100,
            new DepthRange(20, 20), 100
    );

    private static final Map<DepthRange,Integer> MIXED_DEPTH_COUNTS_PART_2 = Map.of(
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
            createDriverGene(SOMATIC_NO_COVERAGE_GENE, true, false),
            createDriverGene(SOMATIC_EXCLUDED_GENE, true, false)
    );

    private static final List<ExonDepth> EXON_DEPTHS = List.of(
            new ExonDepth(GERMLINE_GENE, "1", 100, 149, 1, 0, 0, MIXED_DEPTH_COUNTS_PART_1), // Exons summed per gene
            new ExonDepth(GERMLINE_GENE, "1", 150, 200, 2, 0, 0, MIXED_DEPTH_COUNTS_PART_2),
            new ExonDepth(GERMLINE_CHR_X_GENE, "X", 100, 200, 1, 0, 0, MIXED_DEPTH_COUNTS),
            new ExonDepth(SOMATIC_GENE, "1", 300, 400, 1, 0, 0, MIXED_DEPTH_COUNTS),
            new ExonDepth(SOMATIC_GENE, "1", 450, 460, 2, 0, 0, NO_COVERAGE_DEPTH_COUNTS), // Excluded exon
            new ExonDepth(SOMATIC_CHR_X_GENE, "X", 300, 400, 1, 0, 0, MIXED_DEPTH_COUNTS),
            new ExonDepth(SOMATIC_NO_COVERAGE_GENE, "1", 500, 600, 1, 0, 0, NO_COVERAGE_DEPTH_COUNTS),
            new ExonDepth(SOMATIC_EXCLUDED_GENE, "1", 700, 800, 1, 0, 0, NO_COVERAGE_DEPTH_COUNTS) // All exons excluded
    );

    private static final ExonExclusions EXON_EXCLUSIONS = new ExonExclusions(Map.of(
            SOMATIC_GENE, Set.of(2),
            SOMATIC_EXCLUDED_GENE, Set.of(1)
    ));

    @Test
    public void canGetGermlineMissedVariantLikelihoods()
    {
        List<String> reportableGenes = ReportableGenes.forSampleType(DRIVER_GENES, SampleType.NORMAL);
        List<Feature> actualFeatures = MissedGeneVariantPrep.getMissedVariantLikelihoods(EXON_DEPTHS, reportableGenes, EXON_EXCLUSIONS, 1, Gender.FEMALE);

        List<Feature> expectedFeatures = List.of(
                createExpectedFeature(GERMLINE_GENE, 0.000389385),
                createExpectedFeature(GERMLINE_CHR_X_GENE, 0.000389385) // Diploid chrX in females
        );

        assertFeaturesEqual(expectedFeatures, actualFeatures);
    }

    @Test
    public void canGetSomaticMissedVariantLikelihoods()
    {
        List<String> reportableGenes = ReportableGenes.forSampleType(DRIVER_GENES, SampleType.TUMOR);
        List<Feature> actualFeatures = MissedGeneVariantPrep.getMissedVariantLikelihoods(EXON_DEPTHS, reportableGenes, EXON_EXCLUSIONS, 0.5, Gender.MALE);

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
