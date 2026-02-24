package com.hartwig.hmftools.qsee.prep.category;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.driver.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

import org.junit.Test;

public class MissedGeneVariantPrepTest
{
    private static final String GENE_AR = "AR";
    private static final String GENE_BRCA1 = "BRCA1";

    private static final double DUMMY_MISSED_VARIANT_LIKELIHOOD = 0.0;

    private static final List<DriverGene> DRIVER_GENES = createDriverCatalog();
    private static final List<GeneDepth> GENE_DEPTHS = createGeneDepths();

    @Test
    public void canGetGermlineReportableGenes()
    {
        List<String> reportableGenes = MissedGeneVariantPrep.getReportableGenes(DRIVER_GENES, SampleType.NORMAL);
        List<Feature> actualFeatures = MissedGeneVariantPrep.getMissedVariantLikelihoods(GENE_DEPTHS, reportableGenes);

        List<Feature> expectedFeatures = List.of(
                createExpectedFeature(GENE_BRCA1)
        );

        for(int i = 0; i < expectedFeatures.size(); i++)
        {
            Feature actualFeature = actualFeatures.get(i);
            Feature expectedFeature = expectedFeatures.get(i);
            assertEquals(expectedFeature.key().name(), actualFeature.key().name());
            assertEquals(expectedFeature.value(), actualFeature.value(), 0.001);
        }
    }

    @Test
    public void canGetSomaticReportableGenes()
    {
        List<String> reportableGenes = MissedGeneVariantPrep.getReportableGenes(DRIVER_GENES, SampleType.TUMOR);
        List<Feature> actualFeatures = MissedGeneVariantPrep.getMissedVariantLikelihoods(GENE_DEPTHS, reportableGenes);

        List<Feature> expectedFeatures = List.of(
                createExpectedFeature(GENE_AR),
                createExpectedFeature(GENE_BRCA1)
        );

        for(int i = 0; i < expectedFeatures.size(); i++)
        {
            Feature actualFeature = actualFeatures.get(i);
            Feature expectedFeature = expectedFeatures.get(i);
            assertEquals(expectedFeature.key().name(), actualFeature.key().name());
            assertEquals(expectedFeature.value(), actualFeature.value(), 0.001);
        }
    }

    private static List<DriverGene> createDriverCatalog()
    {
        DriverGene geneAR = ImmutableDriverGene.builder()
                .gene(GENE_AR)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportHetDeletion(false)
                .reportDisruption(false)
                .hetDeletionThreshold(0.0)
                .reportLoh(false)
                .reportAmplification(true)
                .amplificationRatio(0.0)
                .reportSomaticHotspot(true)
                .likelihoodType(DriverCategory.ONCO)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
                .reportGermlineAmplification(false)
                .reportPGX(false)
                .build();

        DriverGene geneBRCA1 = ImmutableDriverGene.builder()
                .gene(GENE_BRCA1)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportHetDeletion(false)
                .reportDisruption(true)
                .hetDeletionThreshold(0.0)
                .reportLoh(false)
                .reportAmplification(false)
                .amplificationRatio(0.0)
                .reportSomaticHotspot(true)
                .likelihoodType(DriverCategory.TSG)
                .reportGermlineVariant(DriverGeneGermlineReporting.ANY)
                .reportGermlineHotspot(DriverGeneGermlineReporting.ANY)
                .reportGermlineDisruption(DriverGeneGermlineReporting.ANY)
                .reportGermlineDeletion(DriverGeneGermlineReporting.ANY)
                .reportGermlineAmplification(false)
                .reportPGX(false)
                .build();

        return List.of(geneAR, geneBRCA1);
    }

    private static List<GeneDepth> createGeneDepths()
    {
        GeneDepth geneAR = new GeneDepth(
                GENE_AR,
                "chrX",
                67545147,
                67723841,
                DUMMY_MISSED_VARIANT_LIKELIHOOD,
                new int[]{ }
        );

        GeneDepth geneBRCA1 = new GeneDepth(
                GENE_BRCA1,
                "chr17",
                43045678,
                43124096,
                DUMMY_MISSED_VARIANT_LIKELIHOOD,
                new int[]{ }
        );

        return List.of(geneAR, geneBRCA1);
    }

    private static Feature createExpectedFeature(String geneName)
    {
        return new Feature("Gene=" + geneName, 0.0, FeatureType.MISSED_VARIANT_LIKELIHOOD, SourceTool.BAM_METRICS, null);
    }
}
