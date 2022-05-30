package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.protect.ProtectEventGenerator;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class Variants {

    private static final Logger LOGGER = LogManager.getLogger(Variants.class);

    private Variants() {
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
        // Reuse PROTECT formatting for ORANGE report.
        String addon = Strings.EMPTY;
        if (!variant.isCanonical()) {
            addon = " (alt)";
        }
        return variant.gene() + addon + " " + AminoAcids.forceSingleLetterProteinAnnotation(ProtectEventGenerator.variantEvent(variant));
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
        if (variant.rnaAlleleReadCount() == null || variant.rnaTotalReadCount() == null) {
            return ReportResources.NOT_AVAILABLE;
        }

        return variant.rnaAlleleReadCount() + " / " + variant.rnaTotalReadCount();
    }
}
