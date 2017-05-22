package com.hartwig.hmftools.common.freec;

import com.hartwig.hmftools.common.copynumber.CopyNumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class FreecCopyNumber implements CopyNumber {

    @Value.Default
    public String genotype() {
        return "-";
    }

    @Value.Default
    public FreecStatus status() {
        return FreecStatus.UNKNOWN;
    }
}
