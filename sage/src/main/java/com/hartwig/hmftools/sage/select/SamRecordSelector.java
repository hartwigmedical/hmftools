package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.utils.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class SamRecordSelector<P extends GenomePosition> extends PositionSelector<P> {

    public SamRecordSelector(@NotNull final List<P> positions) {
        super(positions);
    }

    public void select(final SAMRecord record, final Consumer<P> handler) {
        super.select(record.getAlignmentStart() - SAMRecords.leftSoftClip(record),
                record.getAlignmentEnd() + SAMRecords.rightSoftClip(record),
                handler);
    }
}
