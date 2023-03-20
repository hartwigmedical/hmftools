package com.hartwig.hmftools.datamodel.isofox;

import org.immutables.value.Value;

@Value.Immutable
public interface IsofoxRnaStatistics {
    long totalFragments();

    long duplicateFragments();

    String qcStatus();
}
