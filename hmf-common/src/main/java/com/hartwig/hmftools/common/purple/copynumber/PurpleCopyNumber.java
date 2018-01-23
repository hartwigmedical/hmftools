package com.hartwig.hmftools.common.purple.copynumber;

import com.google.common.base.Strings;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleCopyNumber implements CopyNumber {

    public abstract int bafCount();

    public abstract double averageActualBAF();

    public abstract double averageObservedBAF();

    public abstract double averageTumorCopyNumber();

    public abstract SegmentSupport segmentStartSupport();

    public abstract SegmentSupport segmentEndSupport();

    public abstract CopyNumberMethod method();

    public double minorAllelePloidy() {
        return Doubles.lessThan(averageActualBAF(), 0.50) ? 0 : Math.max(0, (1 - averageActualBAF()) * averageTumorCopyNumber());
    }

    @NotNull
    public String descriptiveBAF() {
        int copyNumber = value();

        int minorAlleleCount = (int) Math.round(minorAllelePloidy());
        int majorAlleleCount = copyNumber - minorAlleleCount;
        return formatBafField("A", Math.max(minorAlleleCount, majorAlleleCount)) + formatBafField("B",
                Math.min(minorAlleleCount, majorAlleleCount));
    }

    @Override
    public int value() {
        return (int) Math.max(0, Math.round(averageTumorCopyNumber()));
    }

    @NotNull
    private static String formatBafField(@NotNull final String allele, final int count) {
        return count < 10 ? Strings.repeat(allele, count) : allele + "[" + count + "x]";
    }
}