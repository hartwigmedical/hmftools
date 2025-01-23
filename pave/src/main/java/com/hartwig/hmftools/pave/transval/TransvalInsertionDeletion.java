package com.hartwig.hmftools.pave.transval;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class TransvalInsertionDeletion extends  TransvalVariant
{
    @NotNull
    final private String DeletedBases;

    @NotNull
    final private Set<TransvalHotspot> Hotspots;
    public TransvalInsertionDeletion(
            @NotNull final String transcriptId,
            @NotNull final String chromosome,
            final int position,
            final boolean spansMultipleExons,
            @NotNull final String referenceNucleotides,
            @NotNull final String deletedBases,
            @NotNull final Set<TransvalHotspot> hotspots)
    {
        super(transcriptId, chromosome, position, spansMultipleExons, referenceNucleotides);
        Hotspots = hotspots;
        Preconditions.checkArgument(referenceNucleotides.contains(deletedBases));
        DeletedBases = deletedBases;
    }

    @NotNull
    public Set<TransvalHotspot> hotspots()
    {
        return Hotspots;
    }

    @NotNull
    public String deletedBases()
    {
        return DeletedBases;
    }

    public int deletedBasesCount()
    {
        return DeletedBases.length();
    }

    @NotNull
    public List<String> alternateCodonsList(final int i)
    {
        return new ArrayList<>();
    }

    public int alternateCodonsCount()
    {
        return 0;
    }
}
