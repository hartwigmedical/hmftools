package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

public class SingleAminoAcidVariant
{
    @NotNull
    final GeneData gene;
    @NotNull
    final TranscriptData transcript;
    @NotNull
    final TranscriptAminoAcids aminoAcidSequence;
    int position;
    @NotNull final String variant;

    private static boolean isValidAminoAcidName(String s)
    {
        // TODO do we need to handle Z and B???
        return AminoAcids.AMINO_ACID_TO_CODON_MAP.containsKey(s);
    }

    public SingleAminoAcidVariant(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position,
            @NotNull final String variant)
    {
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position < transcript.length());
        Preconditions.checkArgument(isValidAminoAcidName(variant));
        this.gene = gene;
        this.transcript = transcript;
        this.aminoAcidSequence = aminoAcidSequence;
        this.position = position;
        this.variant = variant;
    }

    public String referenceAminoAcid()
    {
        return aminoAcidSequence.AminoAcids.substring(position - 1, position);
    }

    public String variantAminoAcid()
    {
        return variant;
    }
}
