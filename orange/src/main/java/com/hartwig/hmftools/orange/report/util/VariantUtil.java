package com.hartwig.hmftools.orange.report.util;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.variant.ReportableVariant;

import org.jetbrains.annotations.NotNull;

public final class VariantUtil {

    private VariantUtil() {
    }

    @NotNull
    public static List<ReportableVariant> sort(@NotNull List<ReportableVariant> variants) {
        return variants.stream().sorted((variant1, variant2) -> {
            if (Math.abs(variant1.driverLikelihood() - variant2.driverLikelihood()) > 0.001) {
                return (variant1.driverLikelihood() - variant2.driverLikelihood()) < 0 ? 1 : -1;
            } else {
                if (variant1.gene().equals(variant2.gene())) {
                    // sort on codon position if gene is the same
                    return extractCodonField(variant1.canonicalHgvsCodingImpact()) - extractCodonField(variant2.canonicalHgvsCodingImpact())
                            < 0 ? -1 : 1;
                } else {
                    return variant1.gene().compareTo(variant2.gene());
                }
            }
        }).collect(Collectors.toList());
    }

    private static int extractCodonField(@NotNull String hgvsCoding) {
        StringBuilder codonAppender = new StringBuilder();
        boolean noDigitFound = true;
        // hgvsCoding starts with "c.", we need to skip that...
        int index = 2;
        while (noDigitFound && index < hgvsCoding.length()) {
            if ((Character.toString(hgvsCoding.charAt(index)).equals("-") && index == 2) || Character.isDigit(hgvsCoding.charAt(index))) {
                codonAppender.append(hgvsCoding.charAt(index));
            } else {
                noDigitFound = false;
            }
            index++;
        }
        return Integer.parseInt(codonAppender.toString());
    }

    @NotNull
    public static String variantField(@NotNull ReportableVariant variant) {
        String consequence = !variant.canonicalHgvsProteinImpact().isEmpty()
                ? AminoAcids.forceSingleLetterProteinAnnotation(variant.canonicalHgvsProteinImpact())
                : variant.canonicalHgvsCodingImpact();
        return variant.gene() + " " + consequence;
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
}
