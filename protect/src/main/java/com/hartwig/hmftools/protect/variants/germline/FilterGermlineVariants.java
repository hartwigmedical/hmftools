package com.hartwig.hmftools.protect.variants.germline;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.protect.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.variants.somatic.DriverSomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class FilterGermlineVariants {

    private FilterGermlineVariants() {
    }

    @NotNull
    public static List<DriverGermlineVariant> filterGermlineVariantsForReporting(@NotNull GermlineReportingModel germlineReportingModel,
            @NotNull List<ReportableGermlineVariant> germlineVariants, @NotNull List<DriverSomaticVariant> driverSomaticVariants,
            @NotNull List<ReportableGainLoss> reportableGainLosses,
            @NotNull List<ReportableHomozygousDisruption> reportableHomozygousDisruptions,
            @NotNull List<ReportableGeneDisruption> reportableGeneDisruptions) {
        Set<String> reportableGermlineGenes = germlineReportingModel.reportableGermlineGenes();

        Set<String> genesWithSomaticDriverMutation =
                driverSomaticVariants.stream().map(x -> x.variant().gene()).collect(Collectors.toSet());
        Set<String> genesWithCopyLoss = reportableGainLosses.stream()
                .filter(reportableGainLoss -> reportableGainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS
                        || reportableGainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS)
                .map(x -> x.gene())
                .collect(Collectors.toSet());
        Set<String> genesWithHomozygousDisruption = reportableHomozygousDisruptions.stream().map(x -> x.gene()).collect(Collectors.toSet());
        Set<String> genesWithGeneDisruption = reportableGeneDisruptions.stream().map(x -> x.gene()).collect(Collectors.toSet());

        Set<String> genesWithSomaticInactivationEvent = Sets.newHashSet();
        genesWithSomaticInactivationEvent.addAll(genesWithSomaticDriverMutation);
        genesWithSomaticInactivationEvent.addAll(genesWithCopyLoss);
        genesWithSomaticInactivationEvent.addAll(genesWithHomozygousDisruption);
        genesWithSomaticInactivationEvent.addAll(genesWithGeneDisruption);

        return filterGermlineVariantsForReporting(reportableGermlineGenes, germlineVariants, genesWithSomaticInactivationEvent);
    }

    @NotNull
    public static List<DriverGermlineVariant> filterGermlineVariantsForReporting(@NotNull Set<String> reportableGermlineGenes,
            @NotNull List<ReportableGermlineVariant> germlineVariants, @NotNull Set<String> genesWithSomaticInactivationEvent) {
        List<DriverGermlineVariant> reportableGermlineVariants = Lists.newArrayList();

        for (ReportableGermlineVariant germlineVariant : presentInTumor(germlineVariants)) {
            if (reportableGermlineGenes.contains(germlineVariant.gene())) {
                // TODO: This will be cleaned up in upcoming germline variant reporting logic
                if (germlineVariant.gene().equals("KIT")) {
                    reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                } else {
                    // Only report germline variants on TSGs if there is a 2nd reportable hit
                    boolean filterBiallelic = germlineVariant.biallelic();

                    boolean filterGermlineVariantInSameGene = false;
                    for (ReportableGermlineVariant variant : germlineVariants) {
                        if (variant != germlineVariant && variant.gene().equals(germlineVariant.gene())) {
                            filterGermlineVariantInSameGene = true;
                        }
                    }

                    boolean filterSomaticVariantInSameGene = genesWithSomaticInactivationEvent.contains(germlineVariant.gene());

                    if (filterBiallelic || filterGermlineVariantInSameGene || filterSomaticVariantInSameGene) {
                        reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                    }
                }
            }
        }

        return reportableGermlineVariants;
    }

    @NotNull
    private static List<ReportableGermlineVariant> presentInTumor(@NotNull List<ReportableGermlineVariant> variants) {
        List<ReportableGermlineVariant> variantsInTumor = Lists.newArrayList();
        for (ReportableGermlineVariant variant : variants) {
            if (isPresentInTumor(variant)) {
                variantsInTumor.add(variant);
            }
        }
        return variantsInTumor;
    }

    private static boolean isPresentInTumor(@NotNull ReportableGermlineVariant germlineVariant) {
        return germlineVariant.adjustedCopyNumber() * germlineVariant.adjustedVaf() >= 0.5;
    }

    @NotNull
    private static DriverGermlineVariant reportableGermlineVariantWithDriverLikelihood(@NotNull ReportableGermlineVariant germlineVariant,
            double driverLikelihood) {
        return ImmutableDriverGermlineVariant.builder().variant(germlineVariant).driverLikelihood(driverLikelihood).build();
    }
}
