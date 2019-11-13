package com.hartwig.hmftools.common.amber;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class BaseDepthFactory {

    private final int minBaseQuality;

    BaseDepthFactory(final int minBaseQuality) {
        this.minBaseQuality = minBaseQuality;
    }

    @NotNull
    public static ModifiableBaseDepth create(@NotNull final BaseDepth pos) {
        return ModifiableBaseDepth.create().from(pos).setIndelCount(0).setRefSupport(0).setAltSupport(0).setReadDepth(0);
    }

    @NotNull
    public static ModifiableBaseDepth create(@NotNull final AmberSite site) {
        return ModifiableBaseDepth.create()
                .setChromosome(site.chromosome())
                .setPosition(site.position())
                .setRef(BaseDepth.Base.valueOf(site.ref()))
                .setAlt(BaseDepth.Base.valueOf(site.alt()))
                .setAltSupport(0)
                .setRefSupport(0)
                .setIndelCount(0)
                .setReadDepth(0);
    }

    void addEvidence(@NotNull final ModifiableBaseDepth evidence, @NotNull final SAMRecord samRecord) {
        int quality = getBaseQuality(evidence, samRecord);
        if (quality >= minBaseQuality) {
            evidence.setReadDepth(evidence.readDepth() + 1);

            int bafPosition = (int) evidence.position();
            int readPosition = samRecord.getReadPositionAtReferencePosition(bafPosition);
            if (readPosition != 0) {
                if (!indel(bafPosition, readPosition, samRecord)) {
                    final char baseChar = samRecord.getReadString().charAt(readPosition - 1);
                    final BaseDepth.Base base = BaseDepth.Base.valueOf(String.valueOf(baseChar).toUpperCase());
                    if (base.equals(evidence.ref())) {
                        evidence.setRefSupport(evidence.refSupport() + 1);
                    } else if (base.equals(evidence.alt())) {
                        evidence.setAltSupport(evidence.altSupport() + 1);
                    }
                } else {
                    evidence.setIndelCount(evidence.indelCount() + 1);
                }
            }
        }
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
