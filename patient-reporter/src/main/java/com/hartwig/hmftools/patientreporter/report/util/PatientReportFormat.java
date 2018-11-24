package com.hartwig.hmftools.patientreporter.report.util;

import static com.google.common.base.Strings.repeat;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PatientReportFormat {

    private PatientReportFormat() {
    }

    @NotNull
    public static String formatNullablePercent(final @Nullable String percentage) {
        String formatTumorPercentage;
        if (percentage != null && !percentage.equals("Not Determined")) {
            formatTumorPercentage = formatPercent(Double.parseDouble(percentage));
        } else if (percentage != null && percentage.equals("Not Determined")) {
            formatTumorPercentage = "Not Determined";
        } else {
            formatTumorPercentage = "N/A";
        }
        return formatTumorPercentage;
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
    public static String correctValueForFitReliability(@NotNull String value, boolean hasReliablePurityFit) {
        return hasReliablePurityFit ? value : "N/A";
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

    @NotNull
    public static String readDepthField(@NotNull AllelicDepth allelicDepth) {
        return allelicDepth.alleleReadCount() + " / " + allelicDepth.totalReadCount();
    }

    @NotNull
    public static String ploidyVafField(double adjustedCopyNumber, double minorAllelePloidy, double adjustedVAF) {
        return descriptiveBAF(adjustedCopyNumber, minorAllelePloidy) + " (" + PatientReportFormat.formatPercent(Math.min(1, adjustedVAF))
                + ")";
    }

    @NotNull
    @VisibleForTesting
    static String descriptiveBAF(double adjustedCopyNumber, double minorAllelePloidy) {
        int totalAlleleCount = (int) Math.max(0, Math.round(adjustedCopyNumber));
        int minorAlleleCount = (int) Math.max(0, Math.round(minorAllelePloidy));
        int majorAlleleCount = Math.max(0, totalAlleleCount - minorAlleleCount);

        return formatBAFField("A", Math.max(minorAlleleCount, majorAlleleCount)) + formatBAFField("B",
                Math.min(minorAlleleCount, majorAlleleCount));
    }

    @NotNull
    private static String formatBAFField(@NotNull String allele, int count) {
        return count < 10 ? repeat(allele, count) : allele + "[" + count + "x]";
    }
}
