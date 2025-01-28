package com.hartwig.hmftools.pave.transval;

import java.util.Set;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class TransvalInsertionDeletion extends  TransvalVariant
{
    @NotNull
    final private String DeletedBases;

    public TransvalInsertionDeletion(
            @NotNull final String transcriptId,
            @NotNull final String chromosome,
            final int position,
            final boolean spansMultipleExons,
            @NotNull final String referenceNucleotides,
            @NotNull final String deletedBases,
            @NotNull final Set<TransvalHotspot> hotspots)
    {
        super(transcriptId, chromosome, position, spansMultipleExons, referenceNucleotides, hotspots);
        Preconditions.checkArgument(referenceNucleotides.contains(deletedBases));
        DeletedBases = deletedBases;
    }

    public int deletedBasesCount()
    {
        return DeletedBases.length();
    }
}
