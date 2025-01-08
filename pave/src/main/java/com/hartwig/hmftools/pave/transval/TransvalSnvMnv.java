package com.hartwig.hmftools.pave.transval;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

import org.jetbrains.annotations.NotNull;

public class TransvalSnvMnv
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
    public final String AlternateNucleotides;
    @NotNull
    public final String ReferenceCodon;
    @NotNull
    public final List<String> AlternateCodons;

    public TransvalSnvMnv(
            @NotNull final String transcriptId,
            @NotNull final String chromosome,
            final int position,
            final boolean spansMultipleExons,
            @NotNull final String referenceNucleotides,
            @NotNull final String alternateNucleotides,
            @NotNull final String referenceCodon,
            @NotNull final List<String> alternateCodons)
    {
        this.TranscriptId = transcriptId;
        Chromosome = RefGenomeFunctions.stripChrPrefix(chromosome);
        Position = position;
        SpansMultipleExons = spansMultipleExons;
        ReferenceNucleotides = referenceNucleotides;
        AlternateNucleotides = alternateNucleotides;
        ReferenceCodon = referenceCodon;
        AlternateCodons = alternateCodons;
    }
}
