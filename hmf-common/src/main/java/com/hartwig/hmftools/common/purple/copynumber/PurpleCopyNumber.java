package com.hartwig.hmftools.common.purple.copynumber;

import com.google.common.base.Strings;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleCopyNumber implements CopyNumber {

    public abstract int bafCount();

    public abstract boolean ratioSupport();

    public abstract double averageActualBAF();

    public abstract double averageObservedBAF();

    public abstract double averageTumorCopyNumber();

    public abstract SegmentSupport support();

    public String descriptiveBAF() {

        int copyNumber = value();
        double constrainedBaf = Math.max(0, Math.min(1, averageActualBAF()));

        int betaAlleleCount = (int) Math.round(constrainedBaf * copyNumber);
        int alphaAlleleCount = copyNumber - betaAlleleCount;
        return formatBafField("A", Math.max(alphaAlleleCount, betaAlleleCount)) + formatBafField("B",
                Math.min(alphaAlleleCount, betaAlleleCount));
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