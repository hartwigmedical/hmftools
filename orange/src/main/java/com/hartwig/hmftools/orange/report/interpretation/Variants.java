package com.hartwig.hmftools.orange.report.interpretation;

import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class Variants {

    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    private Variants() {
    }

    @NotNull
    public static List<VariantEntry> sort(@NotNull List<VariantEntry> variants) {
        return variants.stream().sorted((variant1, variant2) -> {
            double driverLikelihood1 = variant1.driverLikelihood() != null ? variant1.driverLikelihood() : -1;
            double driverLikelihood2 = variant2.driverLikelihood() != null ? variant2.driverLikelihood() : -1;

            int driverCompare = Double.compare(driverLikelihood1, driverLikelihood2);
            if (driverCompare != 0) {
                return driverCompare;
            }

            int geneCompare = variant1.gene().compareTo(variant2.gene());
            if (geneCompare != 0) {
                return geneCompare;
            }

            if (variant1.affectedCodon() == null && variant2.affectedCodon() == null) {
                return 0;
            } else if (variant1.affectedCodon() == null) {
                return 1;
            } else if (variant2.affectedCodon() == null) {
                return -1;
            } else {
                return Integer.compare(variant1.affectedCodon(), variant2.affectedCodon());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static String variantField(@NotNull VariantEntry variant) {
        String addon = Strings.EMPTY;

        if (!variant.isCanonical()) {
            addon = " (alt)";
        }

        return variant.gene() + addon + " " + variant.impact();
    }

    @NotNull
    public static String hotspotField(@NotNull VariantEntry variant) {
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
    public static String rnaDepthField(@NotNull VariantEntry variant) {
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
    public static String driverLikelihoodField(@NotNull VariantEntry variant) {
        return variant.driverLikelihood() != null ? PERCENTAGE_FORMAT.format(variant.driverLikelihood() * 100) : Strings.EMPTY;
    }

    @NotNull
    public static String clonalLikelihoodField(@NotNull VariantEntry variant) {
        return PERCENTAGE_FORMAT.format(100 * variant.clonalLikelihood());
    }

    @NotNull
    public static String phaseSetField(@NotNull VariantEntry variant) {
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
}
