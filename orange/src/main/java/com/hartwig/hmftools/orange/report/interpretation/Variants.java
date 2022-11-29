package com.hartwig.hmftools.orange.report.interpretation;

import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;

import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.orange.algo.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class Variants {

    private static final Logger LOGGER = LogManager.getLogger(Variants.class);

    private static final VariantEffect UPSTREAM_GENE_VARIANT = VariantEffect.UPSTREAM_GENE;

    private static final Set<VariantEffect> PHASED_EFFECTS =
            Sets.newHashSet(VariantEffect.PHASED_INFRAME_DELETION, VariantEffect.PHASED_INFRAME_INSERTION);

    private static final DecimalFormat SINGLE_DIGIT_FORMAT = ReportResources.decimalFormat("#0.0");
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    private Variants() {
    }

    @NotNull
    public static List<PurpleVariant> dedup(@NotNull List<PurpleVariant> variants) {
        List<PurpleVariant> filtered = Lists.newArrayList();
        for (PurpleVariant variant : variants) {
            if (hasCanonicalPhasedEffect(variant) && hasSameEffectWithHigherVCN(variants, variant)) {
                LOGGER.debug("Dedup'ing variant '{}'", variant);
            } else {
                filtered.add(variant);
            }
        }
        return filtered;
    }

    private static boolean hasCanonicalPhasedEffect(@NotNull PurpleVariant variant) {
        for (VariantEffect effect : variant.canonicalImpact().effects()) {
            if (PHASED_EFFECTS.contains(effect)) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasSameEffectWithHigherVCN(@NotNull List<PurpleVariant> variants, @NotNull PurpleVariant variantToMatch) {
        // We assume that variants with same effect have unique hgvs coding impact.
        Double minVariantCopyNumber = null;
        String uniqueHgvsCodingImpact = null;
        PurpleTranscriptImpact variantImpactToMatch = variantToMatch.canonicalImpact();
        for (PurpleVariant variant : variants) {
            PurpleTranscriptImpact variantImpact = variant.canonicalImpact();
            if (variantImpact.effects().equals(variantImpactToMatch.effects()) && variant.gene().equals(variantToMatch.gene())
                    && variantImpact.hgvsProteinImpact().equals(variantImpactToMatch.hgvsProteinImpact())) {
                if (minVariantCopyNumber == null || Doubles.lessThan(variant.variantCopyNumber(), minVariantCopyNumber)) {
                    minVariantCopyNumber = variant.variantCopyNumber();
                    uniqueHgvsCodingImpact = variantImpact.hgvsCodingImpact();
                } else if (Doubles.equal(variant.variantCopyNumber(), minVariantCopyNumber)) {
                    uniqueHgvsCodingImpact = variantImpact.hgvsCodingImpact().compareTo(uniqueHgvsCodingImpact) > 0
                            ? variantImpact.hgvsCodingImpact()
                            : uniqueHgvsCodingImpact;
                }
            }
        }

        boolean matchesMinAlleleCopyNumber = Doubles.equal(variantToMatch.variantCopyNumber(), minVariantCopyNumber);
        boolean matchesBestHgvsCodingImpact = variantImpactToMatch.hgvsCodingImpact().equals(uniqueHgvsCodingImpact);
        return !(matchesMinAlleleCopyNumber && matchesBestHgvsCodingImpact);
    }

    @NotNull
    public static List<PurpleVariant> sort(@NotNull List<PurpleVariant> variants, @NotNull List<DriverCatalog> drivers) {
        return variants.stream().sorted((variant1, variant2) -> {
            DriverCatalog driver1 = Drivers.variantEntryForGene(drivers, variant1.gene());
            DriverCatalog driver2 = Drivers.variantEntryForGene(drivers, variant2.gene());

            double driverLikelihood1 = driver1 != null ? driver1.driverLikelihood() : -1;
            double driverLikelihood2 = driver2 != null ? driver2.driverLikelihood() : -1;

            if (Math.abs(driverLikelihood1 - driverLikelihood2) > 0.001) {
                return (driverLikelihood1 - driverLikelihood2) < 0 ? 1 : -1;
            } else {
                if (variant1.gene().equals(variant2.gene())) {
                    // sort on codon position if gene is the same
                    if (variant1.canonicalImpact().hgvsCodingImpact().isEmpty()) {
                        return 1;
                    } else if (variant2.canonicalImpact().hgvsCodingImpact().isEmpty()) {
                        return -1;
                    } else {
                        int codonVariant1 = extractCodonField(variant1.canonicalImpact().hgvsCodingImpact());
                        int codonVariant2 = extractCodonField(variant2.canonicalImpact().hgvsCodingImpact());
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
    public static String variantField(@NotNull PurpleVariant variant) {
        String addon = Strings.EMPTY;
        // TODO Solve reporting non-canonical variants.
        //        if (!variant.isCanonical()) {
        //            addon = " (alt)";
        //        }

        return variant.gene() + addon + " " + AminoAcids.forceSingleLetterProteinAnnotation(variantEvent(variant));
    }

    @NotNull
    public static String variantCopyNumberField(@NotNull PurpleVariant variant) {
        return SINGLE_DIGIT_FORMAT.format(variant.adjustedCopyNumber() * Math.max(0, Math.min(1, variant.adjustedVAF())));
    }

    @NotNull
    public static String hotspotField(@NotNull PurpleVariant variant) {
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
    public static String rnaDepthField(@NotNull PurpleVariant variant) {
        AllelicDepth rnaDepth = variant.rnaDepth();

        if (rnaDepth == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        String vafAddon = Strings.EMPTY;
        if (rnaDepth.totalReadCount() > 0) {
            double vaf = rnaDepth.alleleReadCount() / (double) rnaDepth.totalReadCount();
            vafAddon = " (" + PERCENTAGE_FORMAT.format(vaf * 100) + ")";
        }

        return rnaDepth.alleleReadCount() + "/" + rnaDepth.totalReadCount() + vafAddon;
    }

    @NotNull
    public static String driverLikelihoodField(@NotNull PurpleVariant variant, @NotNull List<DriverCatalog> drivers) {
        DriverCatalog driver = Drivers.variantEntryForGene(drivers, variant.gene());
        return driver != null ? PERCENTAGE_FORMAT.format(driver.driverLikelihood() * 100) : Strings.EMPTY;
    }

    @NotNull
    public static String clonalLikelihoodField(@NotNull PurpleVariant variant) {
        return PERCENTAGE_FORMAT.format(100 * (1 - variant.subclonalLikelihood()));
    }

    @NotNull
    public static String phaseSetField(@NotNull PurpleVariant variant) {
        List<Integer> localPhaseSets = variant.localPhaseSets();
        if (localPhaseSets == null || localPhaseSets.isEmpty()) {
            return Strings.EMPTY;
        }

        StringJoiner joiner = new StringJoiner(", ");
        for (Integer localPhaseSet : localPhaseSets) {
            joiner.add(String.valueOf(localPhaseSet));
        }
        return joiner.toString();
    }

    @NotNull
    private static String variantEvent(@NotNull PurpleVariant variant) {
        // TODO Solve formatting of reporting of non-canonical events.
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

        if (effects.contains(UPSTREAM_GENE_VARIANT)) {
            return "upstream";
        }

        StringJoiner joiner = new StringJoiner(", ");
        for (VariantEffect effect : effects) {
            joiner.add(effect.effect());
        }
        return joiner.toString();
    }
}
