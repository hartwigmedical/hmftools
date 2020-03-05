package com.hartwig.hmftools.knowledgebasegenerator;

import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AllGenomicEvents {

    @NotNull
    public abstract KnownAmplificationDeletion knownAmplifications();

    @NotNull
    public abstract KnownAmplificationDeletion knownDeletions();
}
