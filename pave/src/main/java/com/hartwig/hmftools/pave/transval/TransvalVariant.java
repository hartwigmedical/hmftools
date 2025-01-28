package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

import org.jetbrains.annotations.NotNull;

public class TransvalVariant
{
    @NotNull
    public final String TranscriptId;
    @NotNull
    public final String Chromosome;
    public final int Position;
    public final boolean SpansMultipleExons;
    @NotNull
    public final String ReferenceNucleotides;
    @NotNull
    protected final Set<TransvalHotspot> Hotspots;

    public TransvalVariant(
            @NotNull final String transcriptId,
            @NotNull final String chromosome,
            final int position,
            final boolean spansMultipleExons,
            @NotNull final String referenceNucleotides, @NotNull final Set<TransvalHotspot> hotspots)
    {
        this.TranscriptId = transcriptId;
        Chromosome = RefGenomeFunctions.stripChrPrefix(chromosome);
        Position = position;
        SpansMultipleExons = spansMultipleExons;
        ReferenceNucleotides = referenceNucleotides;
        Hotspots = hotspots;
    }

    @NotNull
    public Set<TransvalHotspot> hotspots()
    {
        return Hotspots;
    }
}
