package com.hartwig.hmftools.pave.transval;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

import org.jetbrains.annotations.NotNull;

public class TransvalSNV
{
    @NotNull
    public final String TranscriptId;
    @NotNull
    public final String Chromosome;
    public final int Position;
    public final boolean SpansMultipleExons;
    @NotNull
    public final String ReferenceNucleotide;
    @NotNull
    public final String AlternateNucleotide;
    @NotNull
    public final String ReferenceCodon;
    @NotNull
    public final List<String> AlternateCodons;

    public TransvalSNV(
            @NotNull final String transcriptId,
            @NotNull final String chromosome,
            final int position,
            final boolean spansMultipleExons,
            @NotNull final String referenceNucleotide,
            @NotNull final String alternateNucleotide,
            @NotNull final String referenceCodon,
            @NotNull final List<String> alternateCodons)
    {
        this.TranscriptId = transcriptId;
        Chromosome = RefGenomeFunctions.stripChrPrefix(chromosome);
        Position = position;
        SpansMultipleExons = spansMultipleExons;
        ReferenceNucleotide = referenceNucleotide;
        AlternateNucleotide = alternateNucleotide;
        ReferenceCodon = referenceCodon;
        AlternateCodons = alternateCodons;
    }
}
