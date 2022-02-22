package com.hartwig.hmftools.sage.read;

import static java.lang.Math.round;

import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;

import com.hartwig.hmftools.sage.common.RefSequence;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public final class NumberEvents
{
    public static int numberOfEvents(final SAMRecord record, final RefSequence refGenome)
    {
        return (int)round(numberOfEventsRaw(record, refGenome));
    }

    public static double numberOfEventsRaw(final SAMRecord record, final RefSequence refGenome)
    {
        int nm = rawNM(record, refGenome);

        int additionalIndels = 0;
        int softClippedBases = 0;
        for(CigarElement cigarElement : record.getCigar())
        {
            switch(cigarElement.getOperator())
            {
                case D:
                case I:
                    additionalIndels += cigarElement.getLength() - 1;
                    break;
                case S:
                    softClippedBases += cigarElement.getLength();
                    break;
            }
        }

        // 1+ SUM(SoftClip bases) / 10
        double softClipEvents = softClippedBases > 0 ? (1 + softClippedBases / SC_READ_EVENTS_FACTOR) : 0;

        return nm - additionalIndels + softClipEvents;
    }

    public static int rawNM(final SAMRecord record, final RefSequence refGenome)
    {
        Object nm = record.getAttribute("NM");
        if(nm instanceof Integer)
        {
            return (int) nm;
        }

        int offset = refGenome.alignment().Position - refGenome.alignment().Index - 1;
        return SequenceUtil.calculateSamNmTag(record, refGenome.alignment().Bases, offset);
    }

    public static int numberOfEventsWithMNV(double rawNumberEvents, final String ref, final String alt)
    {
        return (int)round(numberOfEventsWithMNVRaw(rawNumberEvents, ref, alt));
    }

    public static double numberOfEventsWithMNVRaw(double rawNumberEvents, final String ref, final String alt)
    {
        // Number of events includes each SNV as an additional event. This unfairly penalises MNVs.
        int differentBases = 0;
        for(int i = 0; i < alt.length(); i++)
        {
            if(alt.charAt(i) != ref.charAt(i))
            {
                differentBases++;
            }
        }

        // We subtract one later when we actually use this value so we need to add one back in here to be consistent with SNVs and INDELs
        return rawNumberEvents - differentBases + 1;
    }
}
