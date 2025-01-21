package com.hartwig.hmftools.pave.transval;

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

    public TransvalVariant(
            @NotNull final String transcriptId,
            @NotNull final String chromosome,
            final int position,
            final boolean spansMultipleExons,
            @NotNull final String referenceNucleotides)
    {
        this.TranscriptId = transcriptId;
        Chromosome = RefGenomeFunctions.stripChrPrefix(chromosome);
        Position = position;
        SpansMultipleExons = spansMultipleExons;
        ReferenceNucleotides = referenceNucleotides;
    }
}
