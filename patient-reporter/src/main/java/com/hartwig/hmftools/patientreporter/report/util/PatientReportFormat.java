package com.hartwig.hmftools.patientreporter.report.util;

import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PatientReportFormat {

    private PatientReportFormat() {
    }

    @NotNull
    public static String formatNullablePercent(final @Nullable Double percentage) {
        return percentage != null ? formatPercent(percentage) : "N/A";
    }

    @NotNull
    public static String formatPercent(final double percentage) {
        return Long.toString(Math.round(percentage * 100D)) + "%";
    }

    @NotNull
    public static String formatPercentWithDefaultCutoffs(double percentage, double minCutoff, double maxCutoff) {
        if (percentage < minCutoff) {
            return "<" + formatPercent(minCutoff);
        } else if (percentage > maxCutoff) {
            return ">" + formatPercent(maxCutoff);
        } else {
            return formatPercent(percentage);
        }
    }

    @NotNull
    public static String ploidyToCopiesString(@Nullable Double ploidy) {
        return ploidy != null ? String.format("%.1f", ploidy) : "-";
    }

    @NotNull
    public static String correctValueForFitStatus(@NotNull FittedPurityStatus fitStatus, @NotNull String value) {
        return fitStatus == FittedPurityStatus.NO_TUMOR ? "N/A" : value;
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
