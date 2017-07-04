package com.hartwig.hmftools.common.purple.copynumber;

import com.google.common.base.Strings;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleCopyNumber implements CopyNumber {

    public static final double CLONAL_DISTANCE = 0.25;

    public abstract int bafCount();

    public abstract boolean ratioSupport();

    public abstract double averageActualBAF();

    public abstract double averageObservedBAF();

    public abstract double averageTumorCopyNumber();

    public abstract StructuralVariantSupport structuralVariantSupport();

    public String descriptiveBAF() {

        int copyNumber = value();
        double constrainedBaf = Math.max(0, Math.min(1, averageActualBAF()));

        int betaAlleleCount = (int) Math.round(constrainedBaf * copyNumber);
        int alphaAlleleCount = copyNumber - betaAlleleCount;
        return Strings.repeat("A", Math.max(alphaAlleleCount, betaAlleleCount)) + Strings.repeat("B",
                Math.min(alphaAlleleCount, betaAlleleCount));
    }

    public boolean isClonal() {
        return Doubles.lessOrEqual(Doubles.distanceFromInteger(averageTumorCopyNumber()), CLONAL_DISTANCE);
    }

    @Override
    public int value() {
        return (int) Math.max(0, Math.round(averageTumorCopyNumber()));
    }
}