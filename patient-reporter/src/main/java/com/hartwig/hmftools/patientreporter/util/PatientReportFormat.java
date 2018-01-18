package com.hartwig.hmftools.patientreporter.util;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum PatientReportFormat {
    ;

    @NotNull
    public static String formatNullablePercent(final @Nullable Double percentage) {
        return percentage == null ? "Na" : formatPercent(percentage);
    }

    @NotNull
    public static String formatPercent(final double percentage) {
        return Long.toString(Math.round(percentage * 100D)) + "%";
    }

    @NotNull
    public static String positionString(@NotNull final GeneAnnotation g) {
        return positionString(g.variant(), g.isStart());
    }

    @NotNull
    public static String positionString(@NotNull final StructuralVariant sv, final boolean start) {
        return String.format("chr%s:%d", sv.chromosome(start), sv.position(start));
    }

    @NotNull
    public static String exonDescription(@NotNull final Transcript t, final boolean upstream) {
        if (t.isPromoter()) {
            return "Promoter Region";
        } else if (t.isExonic()) {
            return String.format("Exon %d", upstream ? t.exonUpstream() : t.exonDownstream());
        } else if (t.isIntronic()) {
            return String.format("Intron %d", t.exonUpstream());
        } else {
            return String.format("Error up(%d) down(%d)", t.exonUpstream(), t.exonDownstream());
        }
    }

    @Nullable
    public static Double alleleFrequency(@NotNull final GeneAnnotation g) {
        return g.isStart() ? g.variant().start().alleleFrequency() : g.variant().end().alleleFrequency();
    }
}
