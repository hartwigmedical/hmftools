package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.utils.sam.SAMRecords;
import com.hartwig.hmftools.sage.sam.SkippedAtLocation;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class SamRecordSelector<P extends GenomePosition> extends PositionSelector<P> {

    private final int maxSkippedReferenceRegions;

    public SamRecordSelector(final int maxSkippedReferenceRegions, @NotNull final List<P> positions) {
        super(positions);
        this.maxSkippedReferenceRegions = maxSkippedReferenceRegions;
    }

    public void select(final SAMRecord record, final Consumer<P> handler) {
        long startWithSoftClip = record.getAlignmentStart() - SAMRecords.leftSoftClip(record);
        long endWithSoftClip = record.getAlignmentEnd() + SAMRecords.rightSoftClip(record);

        final Consumer<P> consumer = position -> {
            if (maxSkippedReferenceRegions > -1 && !SkippedAtLocation.inSkipReference(maxSkippedReferenceRegions, (int) position.position(), record)) {
                handler.accept(position);
            }
        };

        super.select(startWithSoftClip, endWithSoftClip, consumer);
    }

}
