package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.samtools.SamRecordUtils;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class SamRecordSelector<P extends GenomePosition> extends PositionSelector<P>
{
    public SamRecordSelector(@NotNull final List<P> positions)
    {
        super(positions);
    }

    public void select(final SAMRecord record, final Consumer<P> handler)
    {
        int startWithSoftClip = record.getAlignmentStart() - SamRecordUtils.leftSoftClip(record);
        int endWithSoftClip = record.getAlignmentEnd() + SamRecordUtils.rightSoftClip(record);

        super.select(startWithSoftClip, endWithSoftClip, handler);
    }
}
