package com.hartwig.hmftools.sage.context;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.cigar.CigarHandler;
import com.hartwig.hmftools.sage.cigar.CigarTraversal;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ReadContextCounter implements GenomePosition, Consumer<SAMRecord>, CigarHandler {
    private final VariantHotspot hotspot;
    private final ReadContext readContext;

    private int full;
    private int partial;
    private int realigned;

    private int quality;
    private int baseQuality;
    private int mapQuality;

    private int excessiveInferredSize;
    private int inconsistentChromosome;
    private int improperPair;

    private int coverage;

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

    public int[] rcq() {
        return new int[] { improperPair, inconsistentChromosome, excessiveInferredSize };
    }

    public ReadContext readContext() {
        return readContext;
    }

    @Override
    public String toString() {
        return readContext.toString();
    }

    @Override
    public void accept(final SAMRecord record) {
        if (record.getAlignmentStart() <= hotspot.position() && record.getAlignmentEnd() >= hotspot.position()) {
            coverage++;
            if (coverage < 1000) {
                CigarTraversal.traverseCigar(record, this);
            }
        }
    }

    @Override
    public void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readStartIndex,
            final int refStartPosition) {
        byte[] readBases = record.getReadBases();
        for (int i = 0; i < element.getLength(); i++) {
            int readIndex = readStartIndex + i;
            int refPosition = refStartPosition + i;
            if (incrementCounters(refPosition, readIndex, readBases)) {
                incrementQualityScores(readIndex, record);

                if (Math.abs(record.getInferredInsertSize()) >= 1000) {
                    excessiveInferredSize++;
                }

                if (!record.getReferenceName().equals(record.getMateReferenceName())) {
                    inconsistentChromosome++;
                }

                if (!record.getProperPairFlag()) {
                    improperPair++;
                }
            }
        }
    }

    @Override
    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
            final int refPosition) {
        // Empty
    }

    @Override
    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
            final int refPosition) {
        // Empty
    }

    private void incrementQualityScores(int readBasePosition, final SAMRecord record) {
        final int distanceFromReadEdge = Math.min(readBasePosition, record.getReadBases().length - readBasePosition - 1);
        final int mapQuality = record.getMappingQuality();
        final int baseQuality = record.getBaseQualities()[readBasePosition];
        this.mapQuality += mapQuality;
        this.baseQuality += baseQuality;
        this.quality += quality(mapQuality, baseQuality, distanceFromReadEdge);
    }

    private double quality(int mapQuality, int baseQuality, int distanceFromEdge) {
        int quality = Math.min(distanceFromEdge, Math.min(baseQuality - 12, mapQuality - 24));
        return Math.max(0, quality);
    }

    public boolean incrementCounters(long refPosition, int otherReadBytePosition, byte[] otherReadByte) {
        final ReadContext.ReadContextMatch match = readContext.match(otherReadBytePosition, otherReadByte);
        if (!match.equals(ReadContext.ReadContextMatch.NONE)) {
            if (refPosition == hotspot.position()) {
                if (match.equals(ReadContext.ReadContextMatch.FULL)) {
                    full++;
                    return true;
                } else {
                    partial++;
                    return true;
                }
            } else if (match.equals(ReadContext.ReadContextMatch.FULL)) {
                realigned++;
            }
        }
        return false;
    }

}
