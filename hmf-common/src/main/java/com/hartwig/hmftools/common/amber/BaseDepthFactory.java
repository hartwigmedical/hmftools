package com.hartwig.hmftools.common.amber;

import java.util.EnumMap;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

class BaseDepthFactory {

    private final int minBaseQuality;

    BaseDepthFactory(final int minBaseQuality) {
        this.minBaseQuality = minBaseQuality;
    }

    @NotNull
    public static ModifiableBaseDepth create(@NotNull final BaseDepth pos) {
        return ModifiableBaseDepth.create()
                .setChromosome(pos.chromosome())
                .setPosition(pos.position())
                .setRef(pos.ref())
                .setBaseMap(new EnumMap<>(BaseDepth.Base.class))
                .setIndelCount(0)
                .setReadDepth(0);
    }

    @NotNull
    public static ModifiableBaseDepth create(@NotNull final GenomeRegion pos) {
        return ModifiableBaseDepth.create()
                .setChromosome(pos.chromosome())
                .setPosition(pos.start())
                .setBaseMap(new EnumMap<>(BaseDepth.Base.class))
                .setIndelCount(0)
                .setReadDepth(0);
    }

    @NotNull
    public ModifiableBaseDepth addEvidence(@NotNull final ModifiableBaseDepth evidence, @NotNull final SAMRecord samRecord) {
        int quality = getBaseQuality(evidence, samRecord);
        if (quality >= minBaseQuality) {
            evidence.setReadDepth(evidence.readDepth() + 1);

            int bafPosition = (int) evidence.position();
            int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
            if (readPosition != 0) {
                if (!indel(bafPosition, readPosition, samRecord)) {
                    final char baseChar = samRecord.getReadString().charAt(readPosition - 1);
                    final BaseDepth.Base base = BaseDepth.Base.valueOf(String.valueOf(baseChar).toUpperCase());
                    evidence.baseMap().merge(base, 1, (integer, integer2) -> integer + integer2);
                } else {
                    evidence.setIndelCount(evidence.indelCount() + 1);
                }
            }
        }

        return evidence;
    }

    static boolean indel(int bafPosition, int readPosition, @NotNull final SAMRecord samRecord) {

        if (samRecord.getAlignmentEnd() > bafPosition) {

            // Delete?
            if (samRecord.getReadPositionAtReferencePosition(bafPosition + 1) == 0) {
                return true;
            }

            // Insert?
            return samRecord.getReferencePositionAtReadPosition(readPosition + 1) != bafPosition + 1;
        }

        return false;
    }

    static int getBaseQuality(@NotNull final GenomePosition position, @NotNull final SAMRecord samRecord) {
        // Get quality of base after del if necessary
        for (int pos = (int) position.position(); pos <= samRecord.getAlignmentEnd(); pos++) {
            int readPosition = samRecord.getReadPositionAtReferencePosition(pos);
            if (readPosition != 0) {
                return SAMRecords.getBaseQuality(samRecord, readPosition);
            }
        }

        return 0;
    }

}
