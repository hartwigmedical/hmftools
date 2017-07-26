package com.hartwig.hmftools.patientreporter.copynumber;

import com.google.common.annotations.VisibleForTesting;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CopyNumberReport implements Comparable<CopyNumberReport> {

    @VisibleForTesting
    static final String COPY_NUMBER_GAIN = "copy-gain";
    @VisibleForTesting
    static final String COPY_NUMBER_LOSS = "copy-loss";
    @VisibleForTesting
    static final String COPY_NUMBER_NEUTRAL = "none";

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String chromosomeBand();

    @NotNull
    public abstract String gene();

    public abstract int copyNumber();

    @NotNull
    public abstract CopyNumberReportType type();

    @NotNull
    public String description() {
        return type().description();
    }

    //TODO: Once we get rid of freec, this class can implement GenomeRegion and gets comparable for free
    @Override
    public int compareTo(@NotNull final CopyNumberReport other) {
        final Integer intChrom1 = safeToInteger(chromosome());
        final Integer intChrom2 = safeToInteger(other.chromosome());
        if (intChrom1 == null && intChrom2 == null) {
            if (chromosome().equals(other.chromosome())) {
                return chromosomeBand().compareTo(other.chromosomeBand());
            } else {
                return chromosome().compareTo(other.chromosome());
            }
        } else if (intChrom1 == null) {
            return 1;
        } else if (intChrom2 == null) {
            return -1;
        } else {
            if (intChrom1.compareTo(intChrom2) == 0) {
                return chromosomeBand().compareTo(other.chromosomeBand());
            } else {
                return intChrom1.compareTo(intChrom2);
            }
        }
    }

    @Nullable
    private static Integer safeToInteger(@NotNull String string) {
        try {
            return Integer.parseInt(string);
        } catch (NumberFormatException exception) {
            return null;
        }
    }
}
