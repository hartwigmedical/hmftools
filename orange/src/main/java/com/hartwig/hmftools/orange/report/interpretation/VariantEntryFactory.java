package com.hartwig.hmftools.orange.report.interpretation;

import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class VariantEntryFactory {

    private VariantEntryFactory() {
    }

    @NotNull
    public static List<VariantEntry> create(@NotNull List<PurpleVariant> variants, @NotNull List<DriverCatalog> drivers) {
        List<VariantEntry> entries = Lists.newArrayList();
        for (PurpleVariant variant : variants) {
            DriverCatalog driver = Drivers.canonicalMutationEntryForGene(drivers, variant.gene());
            entries.add(toVariantEntry(variant, driver));
        }

        for (DriverCatalog nonCanonicalDriver : Drivers.nonCanonicalMutationEntries(drivers)) {
            PurpleVariant nonCanonicalVariant = findReportedVariantForDriver(variants, nonCanonicalDriver);
            if (nonCanonicalVariant != null) {
                entries.add(toVariantEntry(nonCanonicalVariant, nonCanonicalDriver));
            }
        }

        return entries;
    }

    @NotNull
    private static VariantEntry toVariantEntry(@NotNull PurpleVariant variant, @Nullable DriverCatalog driver) {
        PurpleTranscriptImpact transcriptImpact;

        if (driver != null) {
            transcriptImpact = findTranscriptImpact(variant, driver.transcript());
            if (transcriptImpact == null) {
                throw new IllegalStateException("Could not find impact on transcript " + driver.transcript() + " for variant " + variant);
            }
        } else {
            transcriptImpact = variant.canonicalImpact();
        }

        return ImmutableVariantEntry.builder()
                .gene(variant.gene())
                .isCanonical(driver == null || driver.transcript().equals(variant.canonicalImpact().transcript()))
                .affectedCodon(transcriptImpact.affectedCodon())
                .impact(determineImpact(transcriptImpact))
                .variantCopyNumber(variant.adjustedCopyNumber() * Math.max(0, Math.min(1, variant.adjustedVAF())))
                .totalCopyNumber(variant.adjustedCopyNumber())
                .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                .biallelic(variant.biallelic())
                .hotspot(variant.hotspot())
                .driverLikelihood(driver != null ? driver.driverLikelihood() : null)
                .clonalLikelihood(1 - variant.subclonalLikelihood())
                .localPhaseSets(variant.localPhaseSets())
                .rnaDepth(variant.rnaDepth())
                .genotypeStatus(variant.genotypeStatus())
                .build();
    }

    @Nullable
    private static PurpleVariant findReportedVariantForDriver(@NotNull List<PurpleVariant> variants, @NotNull DriverCatalog driver) {
        List<PurpleVariant> reportedVariantsForGene = findReportedVariantsForGene(variants, driver.gene());
        for (PurpleVariant variant : reportedVariantsForGene) {
            if (findTranscriptImpact(variant, driver.transcript()) != null) {
                return variant;
            }
        }

        return null;
    }

    @NotNull
    private static List<PurpleVariant> findReportedVariantsForGene(@NotNull List<PurpleVariant> variants, @NotNull String geneToFind) {
        List<PurpleVariant> reportedVariantsForGene = Lists.newArrayList();
        for (PurpleVariant variant : variants) {
            if (variant.reported() && variant.gene().equals(geneToFind)) {
                reportedVariantsForGene.add(variant);
            }
        }
        return reportedVariantsForGene;
    }

    @NotNull
    @VisibleForTesting
    static PurpleTranscriptImpact findTranscriptImpact(@NotNull PurpleVariant variant, @NotNull String transcriptToFind) {
        if (variant.canonicalImpact().transcript().equals(transcriptToFind)) {
            return variant.canonicalImpact();
        }

        for (PurpleTranscriptImpact otherImpact : variant.otherImpacts()) {
            if (otherImpact.transcript().equals(transcriptToFind)) {
                return otherImpact;
            }
        }

        return null;
    }

    @NotNull
    @VisibleForTesting
    static String determineImpact(@NotNull PurpleTranscriptImpact impact) {
        String hgvsProteinImpact = impact.hgvsProteinImpact();
        if (!hgvsProteinImpact.isEmpty() && !hgvsProteinImpact.equals("p.?")) {
            return hgvsProteinImpact;
        }

        String hgvsCodingImpact = impact.hgvsCodingImpact();
        if (!hgvsCodingImpact.isEmpty()) {
            return impact.codingEffect() == SPLICE ? hgvsCodingImpact + " splice" : hgvsCodingImpact;
        }

        Set<VariantEffect> effects = impact.effects();
        if (effects.contains(VariantEffect.UPSTREAM_GENE)) {
            return "upstream";
        }

        StringJoiner joiner = new StringJoiner(", ");
        for (VariantEffect effect : effects) {
            joiner.add(effect.effect());
        }
        return joiner.toString();
    }
}
