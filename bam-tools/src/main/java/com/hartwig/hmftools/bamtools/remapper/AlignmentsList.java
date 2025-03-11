package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;
import java.util.stream.Stream;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

import htsjdk.samtools.SAMFlag;

public class AlignmentsList
{
    private final List<BwaMemAlignment> mAlignments;

    public AlignmentsList(final List<BwaMemAlignment> alignments)
    {
        long numberOfNonSupplementaries = alignments.stream().filter(AlignmentsList::isNotSupplementary).count();
        if(numberOfNonSupplementaries != 1)
        {
            throw new IllegalArgumentException("only one non-supplementary alignment is allowed");
        }
        mAlignments = alignments;
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
        return mAlignments.get(0);
    }

    public Stream<HlaAlignment> supplementaryAlignments()
    {
        return mAlignments.subList(1, mAlignments.size()).stream().map(HlaAlignment::new);
    }
}
