package com.hartwig.hmftools.patientreporter.util;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum PatientReportFormat {
    ;

    @NotNull
    public static String formatNullablePercent(final @Nullable Double percentage) {
        return percentage != null ? formatPercent(percentage) : "N/A";
    }

    @NotNull
    public static String formatPercent(final double percentage) {
        return Long.toString(Math.round(percentage * 100D)) + "%";
    }

    @NotNull
    public static String ploidyToCopiesString(@Nullable Double ploidy) {
        return ploidy != null ? String.format("%.1f", ploidy) : "-";
    }

    @NotNull
    public static String correctCopyValueForFitStatus(@NotNull final FittedPurityStatus fitStatus, @NotNull final String copiesValue) {
        return fitStatus == FittedPurityStatus.NO_TUMOR ? "N/A" : copiesValue;
    }

    @NotNull
    public static String exonDescription(@NotNull final Transcript transcript) {
        if (transcript.isPromoter()) {
            return "Promoter Region";
        } else if (transcript.isExonic()) {
            assert transcript.exonUpstream() == transcript.exonDownstream();
            return String.format("Exon %d", transcript.exonUpstream());
        } else if (transcript.isIntronic()) {
            return String.format("Intron %d", transcript.exonUpstream());
        } else {
            return String.format("Error up(%d) down(%d)", transcript.exonUpstream(), transcript.exonDownstream());
        }
    }
}
