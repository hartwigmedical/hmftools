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

//        if (hotspot.alt().equals("T") && hotspot.position() == 144930 && record.getReadName().equals("A00624:8:HHKYHDSXX:4:1571:4508:4758")) {
//            System.out.println("Sdf");
//        }


        if (record.getAlignmentStart() <= hotspot.position() && record.getAlignmentEnd() >= hotspot.position()) {
            coverage++;
            if (coverage < 1000) {

                byte[] readBases = record.getReadBases();
                for (int readBasePosition = 0; readBasePosition < readBases.length; readBasePosition++) {
                    long refPosition = record.getReferencePositionAtReadPosition(readBasePosition + 1);
                    if (incrementCounters(refPosition, readBasePosition, readBases)) {
                        incrementQualityScores(readBasePosition, record);

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
        }
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

//        if (hotspot.position() == refPosition && refPosition == 144930 && hotspot.alt().equals("T")) {
////            System.out.println("sdf");
//            ReadContext myNewReadContext = new ReadContext(otherReadBytePosition, otherReadByte);
//            System.out.println(myNewReadContext);
//
//            if (myNewReadContext.toString().equals("TCAGAGGGAATGCACATCAGTGTGGG")) {
//                System.out.println("sdf");
//            }
//
//
//        }




        ReadContext.ReadContextMatch match = readContext.match(otherReadBytePosition, otherReadByte);
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
