package com.hartwig.hmftools.sage.context;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextCounter implements GenomePosition, Consumer<SAMRecord> {
    private final VariantHotspot hotspot;
    private final ReadContext readContext;

    private int full;
    private int partial;
    private int realigned;

    public ReadContextCounter(@NotNull final VariantHotspot hotspot, @NotNull final ReadContext readContext) {
        assert (readContext.isComplete());
        this.hotspot = hotspot;
        this.readContext = readContext;
    }

    @NotNull
    @Override
    public String chromosome() {
        return hotspot.chromosome();
    }

    @Override
    public long position() {
        return hotspot.position();
    }

    @NotNull
    public ReadContext readContext() {
        return readContext;
    }

    public int full() {
        return full;
    }

    public int partial() {
        return partial;
    }

    public int realigned() {
        return realigned;
    }

    @Override
    public void accept(final SAMRecord record) {
        if (record.getAlignmentStart() <= hotspot.position() && record.getAlignmentEnd() >= hotspot.position()) {

            byte[] readBases = record.getReadBases();
            for (int readBasePosition = 0; readBasePosition < readBases.length; readBasePosition++) {
                long refPosition = record.getReferencePositionAtReadPosition(readBasePosition + 1); //TODO: Check + 1
                final ReadContext refPositionContext = new ReadContext(readBasePosition, readBases);
                accept(refPosition, refPositionContext);
            }
        }
    }

    public boolean accept(long refPosition, ReadContext refPositionContext) {
        ReadContext.ReadContextMatch match = readContext.match(refPositionContext);
        if (!match.equals(ReadContext.ReadContextMatch.NONE)) {
            if (refPosition == hotspot.position()) {
                if (match.equals(ReadContext.ReadContextMatch.FULL)) {
                    full++;
                } else {
                    partial++;
                }
            } else if (match.equals(ReadContext.ReadContextMatch.FULL)) {
                realigned++;
            }
            return true;
        }
        return false;
    }

}
