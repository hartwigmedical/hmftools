package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SamRecordSelector<P extends GenomePosition> extends PositionSelector<P>
{
    public SamRecordSelector(@NotNull final List<P> positions)
    {
        super(positions);
    }

    public void select(final SAMRecord record, final Consumer<P> handler)
    {
        int startWithSoftClip = record.getAlignmentStart() - leftSoftClip(record);
        int endWithSoftClip = record.getAlignmentEnd() + rightSoftClip(record);

        super.select(startWithSoftClip, endWithSoftClip, handler);
    }

    private static int leftSoftClip(@NotNull final SAMRecord record)
    {
        Cigar cigar = record.getCigar();
        if(cigar.numCigarElements() > 0)
        {
            CigarElement firstElement = cigar.getCigarElement(0);
            if(firstElement.getOperator() == CigarOperator.S)
            {
                return firstElement.getLength();
            }
        }

        return 0;
    }

    private static int rightSoftClip(@NotNull final SAMRecord record)
    {
        Cigar cigar = record.getCigar();
        if(cigar.numCigarElements() > 0)
        {
            CigarElement lastElement = cigar.getCigarElement(cigar.numCigarElements() - 1);
            if(lastElement.getOperator() == CigarOperator.S)
            {
                return lastElement.getLength();
            }
        }

        return 0;
    }

}
