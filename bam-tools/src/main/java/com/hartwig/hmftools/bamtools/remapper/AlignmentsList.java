package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;
import java.util.stream.Stream;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFlag;

public class AlignmentsList
{
    @NotNull private final List<BwaMemAlignment> Alignments;

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
        this.Alignments = alignments;
    }

    private static boolean isSupplementary(final BwaMemAlignment alignment)
    {
        return SAMFlag.getFlags(alignment.getSamFlag()).contains(SAMFlag.SUPPLEMENTARY_ALIGNMENT);
    }

    private static boolean isNotSupplementary(final BwaMemAlignment alignment)
    {
        return !isSupplementary(alignment);
    }

    public BwaMemAlignment principalAlignment()
    {
        return Alignments.get(0);
    }

    public Stream<HlaAlignment> supplementaryAlignments()
    {
        return Alignments.subList(1, Alignments.size()).stream().map(HlaAlignment::new);
    }
}
