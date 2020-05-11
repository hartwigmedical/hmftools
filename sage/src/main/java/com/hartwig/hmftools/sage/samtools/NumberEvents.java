package com.hartwig.hmftools.sage.samtools;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class NumberEvents {

    public static int numberOfEvents(@NotNull final SAMRecord record) {
        Object nm = record.getAttribute("NM");
        if (!(nm instanceof Integer)) {
            return 0;
        }

        int additionalIndels = 0;
        for (CigarElement cigarElement : record.getCigar()) {
            switch (cigarElement.getOperator()) {
                case D:
                case I:
                    additionalIndels += (cigarElement.getLength() - 1);
            }
        }

        return (Integer) nm - additionalIndels;
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
