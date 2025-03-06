package com.hartwig.hmftools.pavereverse;

import java.util.Set;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

import org.jetbrains.annotations.NotNull;

public class BaseSequenceVariants
{
    @NotNull
    public final TranscriptData mTranscript;
    @NotNull
    public final String mChromosome;
    @NotNull
    protected final Set<BaseSequenceChange> mChanges;

    public BaseSequenceVariants(
            @NotNull final TranscriptData transcript,
            @NotNull final String chromosome,
            @NotNull final Set<BaseSequenceChange> changes)
    {
        this.mTranscript = transcript;
        mChromosome = RefGenomeFunctions.stripChrPrefix(chromosome);
        mChanges = changes;
    }

    public String transcriptName()
    {
        return mTranscript.TransName;
    }

    @NotNull
    public Set<BaseSequenceChange> changes()
    {
        return mChanges;
    }
}
