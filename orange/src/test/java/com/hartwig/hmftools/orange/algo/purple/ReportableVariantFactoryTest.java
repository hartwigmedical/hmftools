package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.junit.Test;

public class ReportableVariantFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canResolveReportableSomaticVariants() {
        String gene1 = "gene1";
        String gene2 = "gene2";
        SomaticVariant variant1 = SomaticVariantTestFactory.builder().reported(true).gene(gene1).canonicalTranscript("transcript1").build();
        SomaticVariant variant2 =
                SomaticVariantTestFactory.builder().reported(false).gene(gene2).canonicalTranscript("transcript2").build();

        double likelihood = 0.6;
        DriverCatalog driverGene1 =
                DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript1", DriverCategory.ONCO);

        List<ReportableVariant> reportable = ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant1, variant2),
                Lists.newArrayList(driverGene1));

        assertEquals(1, reportable.size());
        assertEquals(likelihood, reportable.get(0).driverLikelihood(), EPSILON);
    }

    @Test
    public void canResolveGermlineVariantsWithMultipleDrivers() {
        String gene = "gene";
        SomaticVariant variant = SomaticVariantTestFactory.builder().reported(true).gene(gene).canonicalTranscript("transcript1").build();

        DriverCatalog driver1 =
                DriverCatalogTestFactory.createCanonicalGermlineMutationEntryForGene(gene, 0.6, "transcript1", DriverCategory.ONCO);
        DriverCatalog driver2 =
                ImmutableDriverCatalog.builder().from(driver1).driver(DriverType.GERMLINE_DELETION).driverLikelihood(1D).build();

        List<ReportableVariant> reportable =
                ReportableVariantFactory.toReportableGermlineVariants(Lists.newArrayList(variant), Lists.newArrayList(driver1, driver2));

        assertEquals(0.6, reportable.get(0).driverLikelihood(), EPSILON);
    }

    @Test
    public void canResolveReportableFromNonCanonicalDrivers() {
        String gene = "gene";
        SomaticVariant variant = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene(gene)
                .canonicalTranscript("transcript1")
                .otherReportedEffects("ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT")
                .build();

        double likelihood = 0.6;
        DriverCatalog driverNonCanonical = DriverCatalogTestFactory.createNonCanonicalSomaticMutationEntryForGene(gene,
                likelihood,
                "ENST00000579755",
                DriverCategory.ONCO);

        List<ReportableVariant> reportable =
                ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant), Lists.newArrayList(driverNonCanonical));

        assertEquals(1, reportable.size());
        assertEquals(likelihood, reportable.get(0).driverLikelihood(), EPSILON);

        double likelihoodCanonical = 0.6;
        SomaticVariant variant2 = SomaticVariantTestFactory.builder()
                .reported(true)
                .gene(gene)
                .canonicalTranscript("transcript2")
                .otherReportedEffects("ENST00000579755|c.246_247delCG|p.Gly83fs|frameshift_variant|NONSENSE_OR_FRAMESHIFT")
                .build();
        DriverCatalog driverCanonical = ImmutableDriverCatalog.builder()
                .from(DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene(gene,
                        likelihoodCanonical,
                        "transcript2",
                        DriverCategory.ONCO))
                .isCanonical(true)
                .build();
        List<ReportableVariant> reportable2 = ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant2),
                Lists.newArrayList(driverNonCanonical, driverCanonical));

        assertEquals(2, reportable2.size());
        assertEquals(likelihoodCanonical, reportable2.get(0).driverLikelihood(), EPSILON);
    }
}