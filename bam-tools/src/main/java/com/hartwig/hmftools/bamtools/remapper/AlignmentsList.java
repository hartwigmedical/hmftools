package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFlag;

public class AlignmentsList
{
    @NotNull public final List<BwaMemAlignment> alignments;

    public AlignmentsList(@NotNull final List<BwaMemAlignment> alignments)
    {
        if (alignments.isEmpty()) {
            throw new IllegalArgumentException("alignments must not be empty");
        }
        if (isSupplementary(alignments.get(0)))
        {
            throw new IllegalArgumentException("first alignment must not be supplementary");
        }
        long numberOfNonSupplementaries = alignments.stream().filter(AlignmentsList::isNotSupplementary).count();
        if (numberOfNonSupplementaries != 1)
        {
            throw new IllegalArgumentException("only one non-supplementary alignment is allowed");
        }
        this.alignments = alignments;
    }

    private static boolean isSupplementary(final BwaMemAlignment alignment)
    {
        return SAMFlag.getFlags(alignment.getSamFlag()).contains(SAMFlag.SUPPLEMENTARY_ALIGNMENT);
    }

    private static boolean isNotSupplementary(final BwaMemAlignment alignment)
    {
        return !isSupplementary(alignment);
    }

    public boolean hasMultipleAlignments()
    {
        return alignments.size() > 1;
    }
}
