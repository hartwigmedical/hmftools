package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.amber.NormalBAFFactory.indel;

import com.hartwig.hmftools.common.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class TumorBAFFactory {

    @NotNull
    public static ModifiableTumorBAF create(@NotNull final NormalBAF normal) {
        return ModifiableTumorBAF.create()
                .from(normal)
                .setRef(normal.ref().toString())
                .setAlt(normal.alt().toString())
                .setNormalReadDepth(normal.readDepth())
                .setNormalRefSupport(normal.refCount())
                .setNormalAltSupport(normal.baseMap().get(normal.alt()))
                .setTumorReadDepth(0)
                .setTumorRefSupport(0)
                .setTumorAltQuality(0)
                .setTumorAltSupport(0);
    }

    @NotNull
    public static ModifiableTumorBAF addEvidence(@NotNull final ModifiableTumorBAF evidence, @NotNull final SAMRecord samRecord) {
        int bafPosition = (int) evidence.position();

        int readPosition = samRecord.getReadPositionAtReferencePosition((int) evidence.position());
        if (readPosition != 0) {
            evidence.setTumorReadDepth(evidence.tumorReadDepth() + 1);
            if (!indel(bafPosition, readPosition, samRecord)) {
                final String base = String.valueOf(samRecord.getReadString().charAt(readPosition - 1));
                if (base.equals(evidence.ref())) {
                    evidence.setTumorRefSupport(evidence.tumorRefSupport() + 1);
                } else if (base.equals(evidence.alt())) {
                    int quality = SAMRecords.getBaseQuality(samRecord.getBaseQualityString().charAt(readPosition - 1));
                    evidence.setTumorAltSupport(evidence.tumorAltSupport() + 1);
                    evidence.setTumorAltQuality(evidence.tumorAltQuality() + quality);
                }
            }
        }
        return evidence;
    }
}
