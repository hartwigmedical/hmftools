package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.copynumber.CopyNumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleCopyNumber implements CopyNumber {

    public abstract int bafCount();

    public abstract double averageObservedBAF();

    public abstract double averageActualBAF();

    public abstract double averageTumorCopyNumber();

    @Override
    public int value() {
        return (int) Math.round(averageTumorCopyNumber());
    }
}