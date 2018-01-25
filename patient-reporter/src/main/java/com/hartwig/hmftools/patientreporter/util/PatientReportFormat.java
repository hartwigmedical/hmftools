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
    public static String positionString(@NotNull final GeneAnnotation geneAnnotation) {
        return positionString(geneAnnotation.variant(), geneAnnotation.isStart());
    }

    @NotNull
    public static String positionString(@NotNull final StructuralVariant sv, final boolean start) {
        return String.format("chr%s:%d", sv.chromosome(start), sv.position(start));
    }

    @NotNull
    public static String exonDescription(@NotNull final Transcript transcript, final boolean upstream) {
        if (transcript.isPromoter()) {
            return "Promoter Region";
        } else if (transcript.isExonic()) {
            return String.format("Exon %d", upstream ? transcript.exonUpstream() : transcript.exonDownstream());
        } else if (transcript.isIntronic()) {
            return String.format("Intron %d", transcript.exonUpstream());
        } else {
            return String.format("Error up(%d) down(%d)", transcript.exonUpstream(), transcript.exonDownstream());
        }
    }
}
