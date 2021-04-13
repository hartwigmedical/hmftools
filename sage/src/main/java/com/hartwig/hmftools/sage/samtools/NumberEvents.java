package com.hartwig.hmftools.sage.samtools;

import com.hartwig.hmftools.sage.ref.RefSequence;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public final class NumberEvents {

    public static int numberOfEvents(@NotNull final SAMRecord record, @NotNull final RefSequence refGenome) {
        int nm = rawNM(record, refGenome);

        int additionalIndels = 0;
        int softClips = 0;
        for (CigarElement cigarElement : record.getCigar()) {
            switch (cigarElement.getOperator()) {
                case D:
                case I:
                    additionalIndels += (cigarElement.getLength() - 1);
                    break;
                case S:
                    softClips++;
                    break;
            }
        }

        return nm - additionalIndels + softClips;
    }

    public static int rawNM(@NotNull final SAMRecord record, @NotNull final RefSequence refGenome) {
        Object nm = record.getAttribute("NM");
        if (nm instanceof Integer) {
            return (int) nm;
        }

        int offset = refGenome.alignment().position() - refGenome.alignment().index() - 1;
        return SequenceUtil.calculateSamNmTag(record, refGenome.alignment().bases(), offset);
    }

    public static int numberOfEventsWithMNV(int rawNumberEvents, @NotNull final String ref, @NotNull final String alt) {
        if (ref.length() == alt.length() && ref.length() > 1) {
            // Number of events includes each SNV as an additional event. This unfairly penalises MNVs.
            int differentBases = 0;
            for (int i = 0; i < alt.length(); i++) {
                if (alt.charAt(i) != ref.charAt(i)) {
                    differentBases++;
                }
            }

            // We subtract one later when we actually use this value so we need to add one back in here to be consistent with SNVs and INDELs
            return rawNumberEvents - differentBases + 1;
        }

        return rawNumberEvents;
    }
}
