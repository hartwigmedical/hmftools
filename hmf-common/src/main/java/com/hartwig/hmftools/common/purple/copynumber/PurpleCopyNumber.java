package com.hartwig.hmftools.common.purple.copynumber;

import com.google.common.base.Strings;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleCopyNumber implements CopyNumber {

    @Value.Default
    public boolean ratioSupport() {
        return true;
    }

    @Value.Default
    public StructuralVariantSupport structuralVariantSupport() {
        return StructuralVariantSupport.NONE;
    }

    public abstract int bafCount();

    public abstract double averageObservedBAF();

    public abstract double averageActualBAF();

    public abstract double averageTumorCopyNumber();

    public String descriptiveBAF() {
        int betaAlleleCount = (int) Math.round(Math.max(0, averageActualBAF()) * value());
        int alphaAlleleCount = value() - betaAlleleCount;
        return Strings.repeat("A", Math.max(alphaAlleleCount, betaAlleleCount)) + Strings.repeat("B",
                Math.min(alphaAlleleCount, betaAlleleCount));
    }

    @Override
    public int value() {
        return (int) Math.max(0, Math.round(averageTumorCopyNumber()));
    }
}