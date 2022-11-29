package com.hartwig.hmftools.orange.report.interpretation;

import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class VariantEntryFactory {

    private static final Logger LOGGER = LogManager.getLogger(VariantEntryFactory.class);

    private VariantEntryFactory() {
    }

    @NotNull
    public static List<VariantEntry> create(@NotNull List<PurpleVariant> variants, @NotNull List<DriverCatalog> drivers) {
        List<VariantEntry> entries = Lists.newArrayList();
        // TODO Add driver likelihood
        for (PurpleVariant variant : variants) {
            entries.add(ImmutableVariantEntry.builder()
                    .gene(variant.gene())
                    .isCanonical(true)
                    .impact(determineImpact(variant))
                    .codon(extractCodonField(variant.canonicalImpact().hgvsCodingImpact()))
                    .variantCopyNumber(variant.adjustedCopyNumber() * Math.max(0, Math.min(1, variant.adjustedVAF())))
                    .totalCopyNumber(variant.adjustedCopyNumber())
                    .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                    .biallelic(variant.biallelic())
                    .hotspot(variant.hotspot())
                    .clonalLikelihood(1 - variant.subclonalLikelihood())
                    .localPhaseSets(variant.localPhaseSets())
                    .rnaDepth(variant.rnaDepth())
                    .genotypeStatus(variant.genotypeStatus())
                    .build());
        }

        // TODO Add non-canonical driver entries.

        return entries;
    }

    @NotNull
    private static String determineImpact(@NotNull PurpleVariant variant) {
        return canonicalVariantEvent(variant);
    }

    @NotNull
    private static String canonicalVariantEvent(@NotNull PurpleVariant variant) {
        return toVariantEvent(variant.canonicalImpact().hgvsProteinImpact(),
                variant.canonicalImpact().hgvsCodingImpact(),
                variant.canonicalImpact().effects(),
                variant.canonicalImpact().codingEffect());
    }

    //    @NotNull
    //    private static String nonCanonicalVariantEvent(@NotNull PurpleVariant variant) {
    //
    //        return toVariantEvent(AltTranscriptReportableInfo.firstOtherHgvsProteinImpact(variant.otherReportedEffects()),
    //                AltTranscriptReportableInfo.firstOtherHgvsCodingImpact(variant.otherReportedEffects()),
    //                AltTranscriptReportableInfo.firstOtherEffects(variant.otherReportedEffects()),
    //                AltTranscriptReportableInfo.firstOtherCodingEffect(variant.otherReportedEffects()));
    //    }

    @NotNull
    @VisibleForTesting
    static String toVariantEvent(@NotNull String protein, @NotNull String coding, @NotNull Set<VariantEffect> effects,
            @NotNull CodingEffect codingEffect) {
        if (!protein.isEmpty() && !protein.equals("p.?")) {
            return protein;
        }

        if (!coding.isEmpty()) {
            return codingEffect == SPLICE ? coding + " splice" : coding;
        }

        if (effects.contains(VariantEffect.UPSTREAM_GENE)) {
            return "upstream";
        }

        StringJoiner joiner = new StringJoiner(", ");
        for (VariantEffect effect : effects) {
            joiner.add(effect.effect());
        }
        return joiner.toString();
    }

    private static int extractCodonField(@NotNull String hgvsCoding) {
        StringBuilder codonAppender = new StringBuilder();
        boolean noDigitFound = true;

        int startIndex = findStartIndex(hgvsCoding);
        int index = startIndex;
        while (noDigitFound && index < hgvsCoding.length()) {
            boolean isMinusSign = Character.toString(hgvsCoding.charAt(index)).equals("-");
            if ((isMinusSign && index == startIndex) || Character.isDigit(hgvsCoding.charAt(index))) {
                codonAppender.append(hgvsCoding.charAt(index));
            } else {
                noDigitFound = false;
            }
            index++;
        }
        String codon = codonAppender.toString();
        if (codon.isEmpty()) {
            LOGGER.warn("Could not extract codon from {}", hgvsCoding);
            return -1;
        } else {
            return Integer.parseInt(codon);
        }
    }

    private static int findStartIndex(@NotNull String hgvsCoding) {
        // hgvsCoding starts with either "c." or "c.*", we need to skip that...
        return hgvsCoding.startsWith("c.*") ? 3 : 2;
    }
}
