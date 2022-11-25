package com.hartwig.hmftools.orange.report.interpretation;

import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;

import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class Variants {

    private static final Logger LOGGER = LogManager.getLogger(Variants.class);

    private static final String UPSTREAM_GENE_VARIANT = "upstream_gene_variant";

    private static final Set<String> PHASED_EFFECTS =
            Sets.newHashSet(VariantEffect.PHASED_INFRAME_DELETION.effect(), VariantEffect.PHASED_INFRAME_INSERTION.effect());

    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    private Variants() {
    }

    @NotNull
    public static List<ReportableVariant> dedup(@NotNull List<ReportableVariant> variants) {
        List<ReportableVariant> filtered = Lists.newArrayList();
        for (ReportableVariant variant : variants) {
            if (PHASED_EFFECTS.contains(variant.canonicalEffect()) && hasSameEffectWithHigherVCN(variants, variant)) {
                LOGGER.debug("Dedup'ing variant '{}'", variant);
            } else {
                filtered.add(variant);
            }
        }
        return filtered;
    }

    private static boolean hasSameEffectWithHigherVCN(@NotNull List<ReportableVariant> variants,
            @NotNull ReportableVariant variantToMatch) {
        // We assume that variants with same effect have unique hgvs coding impact.
        Double minAlleleCopyNumber = null;
        String uniqueHgvsCodingImpact = null;
        for (ReportableVariant variant : variants) {
            if (variant.canonicalEffect().equals(variantToMatch.canonicalEffect()) && variant.gene().equals(variantToMatch.gene())
                    && variant.canonicalHgvsProteinImpact().equals(variantToMatch.canonicalHgvsProteinImpact())) {
                if (minAlleleCopyNumber == null || Doubles.lessThan(variant.alleleCopyNumber(), minAlleleCopyNumber)) {
                    minAlleleCopyNumber = variant.alleleCopyNumber();
                    uniqueHgvsCodingImpact = variant.canonicalHgvsCodingImpact();
                } else if (Doubles.equal(variant.alleleCopyNumber(), minAlleleCopyNumber)) {
                    uniqueHgvsCodingImpact = variant.canonicalHgvsCodingImpact().compareTo(uniqueHgvsCodingImpact) > 0
                            ? variant.canonicalHgvsCodingImpact()
                            : uniqueHgvsCodingImpact;
                }
            }
        }

        boolean matchesMinAlleleCopyNumber = Doubles.equal(variantToMatch.alleleCopyNumber(), minAlleleCopyNumber);
        boolean matchesBestHgvsCodingImpact = variantToMatch.canonicalHgvsCodingImpact().equals(uniqueHgvsCodingImpact);
        return !(matchesMinAlleleCopyNumber && matchesBestHgvsCodingImpact);
    }

    @NotNull
    public static List<ReportableVariant> sort(@NotNull List<ReportableVariant> variants) {
        return variants.stream().sorted((variant1, variant2) -> {
            if (Math.abs(variant1.driverLikelihood() - variant2.driverLikelihood()) > 0.001) {
                return (variant1.driverLikelihood() - variant2.driverLikelihood()) < 0 ? 1 : -1;
            } else {
                if (variant1.gene().equals(variant2.gene())) {
                    // sort on codon position if gene is the same
                    if (variant1.canonicalHgvsCodingImpact().isEmpty()) {
                        return 1;
                    } else if (variant2.canonicalHgvsCodingImpact().isEmpty()) {
                        return -1;
                    } else {
                        int codonVariant1 = extractCodonField(variant1.canonicalHgvsCodingImpact());
                        int codonVariant2 = extractCodonField(variant2.canonicalHgvsCodingImpact());
                        return Integer.compare(codonVariant1, codonVariant2);
                    }
                } else {
                    return variant1.gene().compareTo(variant2.gene());
                }
            }
        }).collect(Collectors.toList());
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

    @NotNull
    public static String variantField(@NotNull ReportableVariant variant) {
        String addon = Strings.EMPTY;
        if (!variant.isCanonical()) {
            addon = " (alt)";
        }

        return variant.gene() + addon + " " + AminoAcids.forceSingleLetterProteinAnnotation(variantEvent(variant));
    }

    @NotNull
    public static String hotspotField(@NotNull ReportableVariant variant) {
        switch (variant.hotspot()) {
            case HOTSPOT:
                return "Yes";
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return "No";
        }
    }

    @NotNull
    public static String rnaDepthField(@NotNull ReportableVariant variant) {
        Integer rnaAlleleReadCount = variant.rnaAlleleReadCount();
        Integer rnaTotalReadCount = variant.rnaTotalReadCount();

        if (rnaAlleleReadCount == null || rnaTotalReadCount == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        String vafAddon = Strings.EMPTY;
        if (rnaTotalReadCount > 0) {
            double vaf = rnaAlleleReadCount / (double) rnaTotalReadCount;
            vafAddon = " (" + PERCENTAGE_FORMAT.format(vaf * 100) + ")";
        }

        return rnaAlleleReadCount + "/" + rnaTotalReadCount + vafAddon;
    }

    @NotNull
    private static String variantEvent(@NotNull ReportableVariant variant) {
        return variant.isCanonical() ? canonicalVariantEvent(variant) : nonCanonicalVariantEvent(variant);
    }

    @NotNull
    private static String canonicalVariantEvent(@NotNull ReportableVariant variant) {
        return toVariantEvent(variant.canonicalHgvsProteinImpact(),
                variant.canonicalHgvsCodingImpact(),
                variant.canonicalEffect(),
                variant.canonicalCodingEffect());
    }

    @NotNull
    private static String nonCanonicalVariantEvent(@NotNull ReportableVariant variant) {
        return toVariantEvent(AltTranscriptReportableInfo.firstOtherHgvsProteinImpact(variant.otherReportedEffects()),
                AltTranscriptReportableInfo.firstOtherHgvsCodingImpact(variant.otherReportedEffects()),
                AltTranscriptReportableInfo.firstOtherEffects(variant.otherReportedEffects()),
                AltTranscriptReportableInfo.firstOtherCodingEffect(variant.otherReportedEffects()));
    }

    @NotNull
    @VisibleForTesting
    static String toVariantEvent(@NotNull String protein, @NotNull String coding, @NotNull String effect,
            @NotNull CodingEffect codingEffect) {
        if (!protein.isEmpty() && !protein.equals("p.?")) {
            return protein;
        }

        if (!coding.isEmpty()) {
            return codingEffect == SPLICE ? coding + " splice" : coding;
        }

        if (effect.equals(UPSTREAM_GENE_VARIANT)) {
            return "upstream";
        }

        return effect;
    }
}
