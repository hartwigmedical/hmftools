package com.hartwig.hmftools.common.cider;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Cdr3Sequence
{
    String cdr3Seq();
    String cdr3AA();
    String locus();
    String filter();
    String blastnStatus();
    int minHighQualBaseReads();
    int assignedReads();
    boolean inFrame();
    boolean containsStop();
}
