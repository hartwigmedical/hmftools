package com.hartwig.hmftools.datamodel.isofox;

import org.immutables.value.Value;

@Value.Immutable
public abstract class IsofoxRnaStatistics {
    public abstract long totalFragments();

    public abstract long duplicateFragments();

    public abstract String qcStatus();
}
