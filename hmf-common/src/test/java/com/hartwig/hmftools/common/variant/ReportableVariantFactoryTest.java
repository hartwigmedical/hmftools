package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogKey;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.test.SomaticVariantTestBuilderFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canMapDriverCatalog() {
        String gene1 = "CDKN2A";
        String gene2 = "BRAF";
        double likelihood = 0.6;

        DriverCatalog driverGene1 = createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript1");
        DriverCatalog driverGene2 = createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript2");
        DriverCatalog driverGene3 = createCanonicalSomaticMutationEntryForGene(gene2, likelihood, "transcript3");
        List<DriverCatalog> mergedDriverCatalog = Lists.newArrayList(driverGene1, driverGene2, driverGene3);


        DriverCatalogKey driverCatalogKey1 = DriverCatalogKey.create(gene1, "transcript1");
        DriverCatalogKey driverCatalogKey2 = DriverCatalogKey.create(gene1, "transcript2");
        DriverCatalogKey driverCatalogKey3 = DriverCatalogKey.create(gene2, "transcript3");

        Map<DriverCatalogKey, DriverCatalog> driverMap = ReportableVariantFactory.toDriverMap(mergedDriverCatalog);

        assertEquals(driverGene1, driverMap.get(driverCatalogKey1));
        assertEquals(driverGene2, driverMap.get(driverCatalogKey2));
        assertEquals(driverGene3, driverMap.get(driverCatalogKey3));
    }

    @Test
    public void canResolveReportableSomaticVariants() {
        String gene1 = "gene1";
        String gene2 = "gene2";
        SomaticVariant variant1 = SomaticVariantTestBuilderFactory.create().reported(true).gene(gene1).build();
        SomaticVariant variant2 = SomaticVariantTestBuilderFactory.create().reported(false).gene(gene2).build();

        double likelihood = 0.6;
        DriverCatalog driverGene1 = createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript1");

        List<ReportableVariant> reportable = ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant1, variant2),
                Lists.newArrayList(driverGene1));

        assertEquals(1, reportable.size());
        assertEquals(likelihood, reportable.get(0).driverLikelihood(), EPSILON);
    }

    @Test
    public void canResolveGermlineVariantsWithMultipleDrivers() {
        String gene = "gene";
        SomaticVariant variant = SomaticVariantTestBuilderFactory.create().reported(true).gene(gene).build();

        DriverCatalog driver1 = createCanonicalGermlineMutationEntryForGene(gene, 0.6, "transcript1");
        DriverCatalog driver2 =
                ImmutableDriverCatalog.builder().from(driver1).driver(DriverType.GERMLINE_DELETION).driverLikelihood(1D).build();

        List<ReportableVariant> reportable =
                ReportableVariantFactory.toReportableGermlineVariants(Lists.newArrayList(variant), Lists.newArrayList(driver1, driver2));

        assertEquals(0.6, reportable.get(0).driverLikelihood(), EPSILON);
    }

    @Test
    public void canResolveReportableFromNonCanonicalDrivers() {
        String gene = "gene";
        SomaticVariant variant = SomaticVariantTestBuilderFactory.create().reported(true).gene(gene).build();

        double likelihood = 0.6;
        DriverCatalog driverNonCanonical = ImmutableDriverCatalog.builder()
                .from(createCanonicalSomaticMutationEntryForGene(gene, likelihood, "transcript1"))
                .isCanonical(false)
                .build();

        List<ReportableVariant> reportable =
                ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant), Lists.newArrayList(driverNonCanonical));

        assertEquals(1, reportable.size());
        assertEquals(likelihood, reportable.get(0).driverLikelihood(), EPSILON);

        double likelihoodCanonical = 0.5;
        SomaticVariant variant2 = SomaticVariantTestBuilderFactory.create().reported(true).gene(gene).build();
        DriverCatalog driverCanonical = createCanonicalSomaticMutationEntryForGene(gene, likelihoodCanonical, "transcript2");
        List<ReportableVariant> reportable2 = ReportableVariantFactory.toReportableSomaticVariants(Lists.newArrayList(variant2),
                Lists.newArrayList(driverNonCanonical, driverCanonical));

        assertEquals(1, reportable2.size());
        assertEquals(likelihoodCanonical, reportable2.get(0).driverLikelihood(), EPSILON);
    }

    @NotNull
    private static DriverCatalog createCanonicalSomaticMutationEntryForGene(@NotNull String gene, double likelihood,
            @NotNull String transcript) {
        return create(gene, likelihood, DriverType.MUTATION, transcript);
    }

    @NotNull
    private static DriverCatalog createCanonicalGermlineMutationEntryForGene(@NotNull String gene, double likelihood,
            @NotNull String transcript) {
        return create(gene, likelihood, DriverType.GERMLINE_MUTATION, transcript);
    }

    private static DriverCatalog create(@NotNull String gene, double likelihood, @NotNull DriverType type, @NotNull String transcript) {
        return ImmutableDriverCatalog.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript(transcript)
                .isCanonical(true)
                .driver(type)
                .category(DriverCategory.ONCO)
                .likelihoodMethod(LikelihoodMethod.DNDS)
                .driverLikelihood(likelihood)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .build();
    }
}
