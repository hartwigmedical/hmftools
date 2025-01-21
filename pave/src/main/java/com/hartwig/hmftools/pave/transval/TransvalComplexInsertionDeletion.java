package com.hartwig.hmftools.pave.transval;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class TransvalComplexInsertionDeletion extends  TransvalVariant
{

    public TransvalComplexInsertionDeletion(
            @NotNull final String transcriptId,
            @NotNull final String chromosome,
            final int position,
            final boolean spansMultipleExons,
            @NotNull final String referenceNucleotides)
    {
        super(transcriptId, chromosome, position, spansMultipleExons, referenceNucleotides);
    }

    public int deletedBasesCount()
    {
        return 0;
    }

    public List<String> alternateCodonsList(final int i)
    {
        return null;
    }

    public int alternateCodonsCount()
    {
        return 0;
    }
}
