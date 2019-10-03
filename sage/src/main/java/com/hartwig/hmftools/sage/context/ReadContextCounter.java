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

    private int quality;
    private int baseQuality;
    private int mapQuality;

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

    public int full() {
        return full;
    }

    public int partial() {
        return partial;
    }

    public int realigned() {
        return realigned;
    }

    public int quality() {
        return quality;
    }

    public int baseQuality() {
        return baseQuality;
    }

    public int mapQuality() {
        return mapQuality;
    }

    public ReadContext readContext() {
        return readContext;
    }


    public void reset() {
        full = 0;
        partial = 0;
        realigned = 0;
        quality = 0;
        baseQuality = 0;
        mapQuality = 0;
    }

    @Override
    public String toString() {
        return readContext.toString();
    }

    @Override
    public void accept(final SAMRecord record) {
        if (record.getAlignmentStart() <= hotspot.position() && record.getAlignmentEnd() >= hotspot.position()) {

            byte[] readBases = record.getReadBases();
            for (int readBasePosition = 0; readBasePosition < readBases.length; readBasePosition++) {
                long refPosition = record.getReferencePositionAtReadPosition(readBasePosition + 1);
                if (incrementCounters(refPosition, readBasePosition, readBases)) {
                    incrementQualityScores(readBasePosition, record);
                }
            }
        }
    }

    private void incrementQualityScores(int readBasePosition, final SAMRecord record) {
        final int distanceFromReadEdge = Math.min(readBasePosition, record.getReadBases().length - readBasePosition - 1);
        final int mapQuality = record.getMappingQuality();
        final int baseQuality = record.getBaseQualities()[readBasePosition];
        final int quality = Math.min(Math.min(mapQuality - 12, baseQuality), distanceFromReadEdge);

        this.mapQuality += mapQuality;
        this.baseQuality += baseQuality;
        this.quality += quality;

    }

    public boolean incrementCounters(long refPosition, int otherReadBytePosition, byte[] otherReadByte) {
        ReadContext.ReadContextMatch match = readContext.match(otherReadBytePosition, otherReadByte);
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
