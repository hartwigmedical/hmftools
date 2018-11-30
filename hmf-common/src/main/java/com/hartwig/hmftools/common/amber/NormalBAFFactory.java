package com.hartwig.hmftools.common.amber;

import java.util.EnumMap;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class NormalBAFFactory {

    private final int minBaseQuality;

    public NormalBAFFactory(final int minBaseQuality) {
        this.minBaseQuality = minBaseQuality;
    }

    @NotNull
    public static ModifiableNormalBAF create(@NotNull final GenomeRegion pos) {
        return ModifiableNormalBAF.create()
                .setChromosome(pos.chromosome())
                .setPosition(pos.start())
                .setBaseMap(new EnumMap<>(NormalBAF.Base.class))
                .setReadDepth(0);
    }

    @NotNull
    public ModifiableNormalBAF addEvidence(@NotNull final ModifiableNormalBAF evidence, @NotNull final SAMRecord samRecord) {
        int bafPosition = (int) evidence.position();

        int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
        if (readPosition != 0) {
            int quality = SAMRecords.getBaseQuality(samRecord, readPosition);
            if (quality >= minBaseQuality) {
                evidence.setReadDepth(evidence.readDepth() + 1);
                if (!indel(bafPosition, readPosition, samRecord)) {
                    final char baseChar = samRecord.getReadString().charAt(readPosition - 1);
                    final NormalBAF.Base base = NormalBAF.Base.valueOf(String.valueOf(baseChar).toUpperCase());
                    evidence.baseMap().merge(base, 1, (integer, integer2) -> integer + integer2);
                }
            }
        }

        return evidence;
    }

    static boolean indel(int bafPosition, int readPosition, @NotNull final SAMRecord samRecord) {

        if (samRecord.getAlignmentEnd() > bafPosition) {

            // Insert?
            if (samRecord.getReadPositionAtReferencePosition(bafPosition + 1) == 0) {
                return true;
            }

            // Delete?
            return samRecord.getReferencePositionAtReadPosition(readPosition + 1) != bafPosition + 1;
        }

        return false;
    }

}
