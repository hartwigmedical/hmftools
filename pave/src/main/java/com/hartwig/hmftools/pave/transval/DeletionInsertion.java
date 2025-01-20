package com.hartwig.hmftools.pave.transval;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

public class DeletionInsertion extends ProteinVariant
{
    final int DeletionLength;
    @NotNull
    final String Alt;

    public DeletionInsertion(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position,
            final int deletionLength,
            @NotNull final String variant)
    {
        super(gene, transcript, aminoAcidSequence, position);
        Preconditions.checkArgument(!variant.isBlank());
        for (char ch : variant.toCharArray())
        {
            Preconditions.checkArgument(isValidAminoAcidName("" + ch));
        }
        this.Position = position;
        this.Alt = variant;
        this.DeletionLength = deletionLength;
    }

    public String altAminoAcidSequence()
    {
        return Alt;
    }

    @Override
    int changedReferenceSequenceLength()
    {
        return DeletionLength;
    }
}
