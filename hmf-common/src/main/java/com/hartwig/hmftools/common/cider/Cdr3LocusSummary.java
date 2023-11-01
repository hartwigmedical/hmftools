package com.hartwig.hmftools.common.cider;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Cdr3LocusSummary {
    String locus();
    int readsUsed();
    int readsTotal();
    boolean downSampled();
    int sequences();
    int passSequences();
}
