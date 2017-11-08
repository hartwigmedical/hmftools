package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.strelka.MNVDetectorApplication.filterMnvRegion;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.strelka.mnv.PotentialMNVRegion;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

final class MNVDetector {

    private MNVDetector() {
    }

    @NotNull
    static List<SomaticVariant> locatePotentialMNVs(@NotNull final List<SomaticVariant> allVariants,
            @NotNull final List<SomaticVariant> reportedVariants) {
        final Set<SomaticVariant> reportedVariantsInMnvs = Sets.newHashSet();
        Pair<PotentialMNVRegion, Optional<PotentialMNVRegion>> outputPair = ImmutablePair.of(PotentialMNVRegion.empty(), Optional.empty());
        final List<VariantContext> variants = allVariants.stream().map(MNVDetector::fromSomaticVariant).collect(Collectors.toList());
        for (final VariantContext variant : variants) {
            final PotentialMNVRegion potentialMNVregion = outputPair.getLeft();
            outputPair = com.hartwig.hmftools.strelka.mnv.MNVDetector.fitsMNVRegion(potentialMNVregion, variant);
            outputPair.getRight()
                    .ifPresent(region -> filterMnvRegion(region).ifPresent(
                            mnvRegion -> reportedVariantsInMnvs.addAll(containedReportedVariants(mnvRegion, reportedVariants))));
        }
        filterMnvRegion(outputPair.getLeft()).ifPresent(
                mnvRegion -> reportedVariantsInMnvs.addAll(containedReportedVariants(mnvRegion, reportedVariants)));
        return reportedVariantsInMnvs.stream().sorted().collect(Collectors.toList());
    }

    @NotNull
    private static Set<SomaticVariant> containedReportedVariants(@NotNull final PotentialMNVRegion mnvRegion,
            @NotNull final List<SomaticVariant> reportedVariants) {
        return reportedVariants.stream()
                .filter(variant -> mnvRegion.variants().stream().anyMatch(mnvVariant -> variantsEqual(mnvVariant, variant)))
                .collect(Collectors.toSet());
    }

    private static boolean variantsEqual(@NotNull final VariantContext mnvVariant, @NotNull final SomaticVariant variant) {
        return mnvVariant.getContig().equals(variant.chromosome()) && mnvVariant.getStart() == variant.position()
                && mnvVariant.getReference().getBaseString().equals(variant.ref());
    }

    @NotNull
    private static VariantContext fromSomaticVariant(@NotNull final SomaticVariant variant) {
        final List<Allele> alleles = buildAlleles(variant);
        return new VariantContextBuilder().chr(variant.chromosome())
                .start(variant.position())
                .stop(variant.position() + variant.ref().length() - 1)
                .alleles(alleles)
                .genotypes(mockGenotype(alleles))
                .make();
    }

    @NotNull
    private static List<Allele> buildAlleles(@NotNull final SomaticVariant variant) {
        final List<Allele> alts = Arrays.stream(variant.alt().split(","))
                .map(String::trim)
                .map(alleleString -> Allele.create(alleleString, false))
                .collect(Collectors.toList());
        final List<Allele> alleles = Lists.newArrayList(Allele.create(variant.ref(), true));
        alleles.addAll(alts);
        return alleles;
    }

    @NotNull
    private static Genotype mockGenotype(@NotNull final List<Allele> alleles) {
        int[] ads = new int[alleles.size()];
        for (int index = 0; index < ads.length; index++) {
            ads[index] = index;
        }
        return new GenotypeBuilder("mockSample").DP(50).AD(ads).make();
    }
}
