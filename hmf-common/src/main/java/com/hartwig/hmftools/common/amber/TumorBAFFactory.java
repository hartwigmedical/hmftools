package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.amber.NormalBAFFactory.indel;

import com.hartwig.hmftools.common.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class TumorBAFFactory {

    private final int minBaseQuality;

    public TumorBAFFactory(final int minBaseQuality) {
        this.minBaseQuality = minBaseQuality;
    }

    @NotNull
    public static ModifiableTumorBAF create(@NotNull final NormalBAF normal) {
        return ModifiableTumorBAF.create()
                .from(normal)
                .setRef(normal.ref().toString())
                .setAlt(normal.alt().toString())
                .setNormalReadDepth(normal.readDepth())
                .setNormalRefSupport(normal.refSupport())
                .setNormalAltSupport(normal.baseMap().get(normal.alt()))
                .setTumorReadDepth(0)
                .setTumorRefSupport(0)
                .setTumorAltQuality(0)
                .setTumorAltSupport(0);
    }

    @NotNull
    public ModifiableTumorBAF addEvidence(@NotNull final ModifiableTumorBAF evidence, @NotNull final SAMRecord samRecord) {
        int bafPosition = (int) evidence.position();

        int readPosition = samRecord.getReadPositionAtReferencePosition((int) evidence.position());
        if (readPosition != 0) {
            int quality = SAMRecords.getBaseQuality(samRecord, readPosition);
            if (quality >= minBaseQuality) {

                evidence.setTumorReadDepth(evidence.tumorReadDepth() + 1);
                if (!indel(bafPosition, readPosition, samRecord)) {
                    final String base = String.valueOf(samRecord.getReadString().charAt(readPosition - 1));
                    if (base.equals(evidence.ref())) {
                        evidence.setTumorRefSupport(evidence.tumorRefSupport() + 1);
                    } else if (base.equals(evidence.alt())) {
                        evidence.setTumorAltSupport(evidence.tumorAltSupport() + 1);
                        evidence.setTumorAltQuality(evidence.tumorAltQuality() + quality);
                    }
                }
            }
        }
        return evidence;
    }
}
