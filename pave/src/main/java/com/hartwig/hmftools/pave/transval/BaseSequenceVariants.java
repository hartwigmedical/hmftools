package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

import org.jetbrains.annotations.NotNull;

public class BaseSequenceVariants
{
    @NotNull
    public final TranscriptData Transcript;
    @NotNull
    public final String Chromosome;
    @NotNull
    protected final Set<BaseSequenceChange> Hotspots;

    public BaseSequenceVariants(
            @NotNull final TranscriptData transcript,
            @NotNull final String chromosome,
            @NotNull final Set<BaseSequenceChange> hotspots)
    {
        this.Transcript = transcript;
        Chromosome = RefGenomeFunctions.stripChrPrefix(chromosome);
        Hotspots = hotspots;
    }

    public String transcriptId()
    {
        return Transcript.TransName;
    }

    @NotNull
    public Set<BaseSequenceChange> hotspots()
    {
        return Hotspots;
    }
}
