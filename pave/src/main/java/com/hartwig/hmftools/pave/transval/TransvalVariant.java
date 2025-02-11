package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

import org.jetbrains.annotations.NotNull;

public class TransvalVariant
{
    @NotNull
    public final TranscriptData Transcript;
    @NotNull
    public final String Chromosome;
//    public final int Position;
    public final boolean SpansMultipleExons;
//    @NotNull
//    public final String ReferenceNucleotides;
    @NotNull
    protected final Set<TransvalHotspot> Hotspots;

    public TransvalVariant(
            @NotNull final TranscriptData transcript,
            @NotNull final String chromosome,
            final boolean spansMultipleExons,
            @NotNull final Set<TransvalHotspot> hotspots)
    {
        this.Transcript = transcript;
        Chromosome = RefGenomeFunctions.stripChrPrefix(chromosome);
//        Position = position;
        SpansMultipleExons = spansMultipleExons;
//        ReferenceNucleotides = referenceNucleotides;
        Hotspots = hotspots;
    }

    public String transcriptId()
    {
        return Transcript.TransName;
    }

    @NotNull
    public Set<TransvalHotspot> hotspots()
    {
        return Hotspots;
    }
}
