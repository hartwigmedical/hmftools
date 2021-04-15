package com.hartwig.hmftools.protect.purple;

import static com.hartwig.hmftools.protect.purple.ReportableVariantFactory.reportableSomaticVariants;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.protect.germline.GermlineReportingEntry;
import com.hartwig.hmftools.protect.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.germline.ImmutableGermlineReportingEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canResolveReportableSomaticVariants() {
        String gene1 = "gene1";
        String gene2 = "gene2";
        SomaticVariant variant1 = SomaticVariantTestBuilderFactory.create().reported(true).gene(gene1).build();
        SomaticVariant variant2 = SomaticVariantTestBuilderFactory.create().reported(false).gene(gene2).build();

        double likelihood = 0.6;
        DriverCatalog driverGene1 = createMutationEntryForGene(gene1, likelihood);
        String notifyGene = "Notify";
        String reportGene = "Report";

        GermlineReportingEntry germlineReportingTrue = ImmutableGermlineReportingEntry.builder()
                .gene(notifyGene)
                .notifyClinicalGeneticist(true)
                .exclusiveHgvsProteinFilter(null)
                .build();

        GermlineReportingEntry germlineReportingFalse = ImmutableGermlineReportingEntry.builder()
                .gene(reportGene)
                .notifyClinicalGeneticist(false)
                .exclusiveHgvsProteinFilter(null)
                .build();
        GermlineReportingModel victim = new GermlineReportingModel(Lists.newArrayList(germlineReportingTrue, germlineReportingFalse));

        List<ReportableVariant> reportable =
                reportableSomaticVariants(Lists.newArrayList(variant1, variant2), Lists.newArrayList(driverGene1), victim);

        assertEquals(1, reportable.size());
        assertEquals(likelihood, reportable.get(0).driverLikelihood(), EPSILON);
    }

    @NotNull
    private static DriverCatalog createMutationEntryForGene(@NotNull String gene, double likelihood) {
        return ImmutableDriverCatalog.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .driver(DriverType.MUTATION)
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
