package com.hartwig.hmftools.patientreporter.report.util;

import static com.google.common.base.Strings.repeat;

import com.hartwig.hmftools.common.variant.AllelicDepth;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PatientReportFormat {

    private PatientReportFormat() {
    }

    @NotNull
    public static String formatPercent(final double percentage) {
        return Long.toString(Math.round(percentage * 100D)) + "%";
    }

    @NotNull
    public static String ploidyToCopiesString(@Nullable Double ploidy) {
        return ploidy != null ? String.format("%.1f", ploidy) : Strings.EMPTY;
    }

    @NotNull
    public static String correctValueForFitReliability(@NotNull String value, boolean hasReliablePurityFit) {
        return hasReliablePurityFit ? value : "N/A";
    }

    @NotNull
    public static String exonDescription(int exonUp, int exonDown) {
        if (exonUp > 0) {
            if (exonUp == exonDown) {
                return String.format("Exon %d", exonUp);
            } else if (exonDown - exonUp == 1) {
                return String.format("Intron %d", exonUp);
            }
        } else if (exonUp == 0 && (exonDown == 1 || exonDown == 2)) {
            return "Promoter Region";
        }

        return String.format("ERROR up=%d, down=%d", exonUp, exonDown);
    }

    @NotNull
    public static String readDepthField(@NotNull AllelicDepth allelicDepth) {
        return allelicDepth.alleleReadCount() + " / " + allelicDepth.totalReadCount();
    }

    @NotNull
    public static String ploidyVafField(double adjustedCopyNumber, double minorAllelePloidy, double adjustedVAF) {
        return descriptiveBAF(adjustedCopyNumber, minorAllelePloidy) + " (" + PatientReportFormat.formatPercent(Math.min(1,
                Math.max(0, adjustedVAF))) + ")";
    }

    @NotNull
    private static String descriptiveBAF(double adjustedCopyNumber, double minorAllelePloidy) {
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
