package com.hartwig.hmftools.protect.variants.germline;

import java.util.List;
import java.util.Map;
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

@Deprecated
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
        Set<String> monoallelicGenesReportable = germlineReportingModel.monoallelicGenesReportable();
        Map<String, String> reportableSpecificVariants = germlineReportingModel.reportableSpecificVariants();

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

        return filterGermlineVariantsForReporting(reportableGermlineGenes,
                germlineVariants,
                genesWithSomaticInactivationEvent,
                monoallelicGenesReportable,
                reportableSpecificVariants);
    }

    @NotNull
    private static List<DriverGermlineVariant> determineReportableVariants(@NotNull Map<String, String> reportableSpecificVariants,
            @NotNull ReportableGermlineVariant germlineVariant, @NotNull List<DriverGermlineVariant> reportableGermlineVariants) {
        boolean containsInReportableSpecificVariants = false;
        for (Map.Entry<String, String> entry : reportableSpecificVariants.entrySet()) {
            if (entry.getValue().contains(germlineVariant.gene())) {
                if (entry.getKey().equals(germlineVariant.hgvsProtein())) {
                    reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
                    containsInReportableSpecificVariants = true;

                }
            }
        }

        if (!containsInReportableSpecificVariants) {
            reportableGermlineVariants.add(reportableGermlineVariantWithDriverLikelihood(germlineVariant, 1.0));
        }
        return reportableGermlineVariants;
    }

    @NotNull
    private static List<DriverGermlineVariant> filterGermlineVariantsForReporting(@NotNull Set<String> reportableGermlineGenes,
            @NotNull List<ReportableGermlineVariant> germlineVariants, @NotNull Set<String> genesWithSomaticInactivationEvent,
            @NotNull Set<String> monoallelicGenesReportable, @NotNull Map<String, String> reportableSpecificVariants) {
        List<DriverGermlineVariant> reportableGermlineVariants = Lists.newArrayList();

        for (ReportableGermlineVariant germlineVariant : presentInTumor(germlineVariants)) {
            if (reportableGermlineGenes.contains(germlineVariant.gene())) {
                if (monoallelicGenesReportable.contains(germlineVariant.gene())) {
                    reportableGermlineVariants =
                            determineReportableVariants(reportableSpecificVariants, germlineVariant, reportableGermlineVariants);
                } else {
                    // Only report germline variants on genes with BIALLELIC condition when a 2nd reportable hit is present
                    boolean filterBiallelic = germlineVariant.biallelic();

                    boolean filterGermlineVariantInSameGene = false;
                    for (ReportableGermlineVariant variant : germlineVariants) {
                        if (variant != germlineVariant && variant.gene().equals(germlineVariant.gene())) {
                            filterGermlineVariantInSameGene = true;
                        }
                    }

                    boolean filterSomaticVariantInSameGene = genesWithSomaticInactivationEvent.contains(germlineVariant.gene());

                    if (filterBiallelic || filterGermlineVariantInSameGene || filterSomaticVariantInSameGene) {
                        reportableGermlineVariants =
                                determineReportableVariants(reportableSpecificVariants, germlineVariant, reportableGermlineVariants);
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

    public static boolean isPresentInTumor(@NotNull ReportableGermlineVariant germlineVariant) {
        return germlineVariant.adjustedCopyNumber() * germlineVariant.adjustedVaf() >= 0.5;
    }

    @NotNull
    private static DriverGermlineVariant reportableGermlineVariantWithDriverLikelihood(@NotNull ReportableGermlineVariant germlineVariant,
            double driverLikelihood) {
        return ImmutableDriverGermlineVariant.builder().variant(germlineVariant).driverLikelihood(driverLikelihood).build();
    }
}
