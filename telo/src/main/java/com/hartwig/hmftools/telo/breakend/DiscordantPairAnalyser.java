package com.hartwig.hmftools.telo.breakend;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.telo.ReadGroup;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

// find all the disordant pairs given
// a list of candidate telomeric break ends
public class DiscordantPairAnalyser
{
    private final int MAX_READ_DISTANCE = 1000;
    private final List<TelomericBreakEnd> mBreakEnds;

    public DiscordantPairAnalyser(List<TelomericBreakEnd> breakEnds)
    {
        mBreakEnds = breakEnds;
    }

    // after we found all the candidate break ends, we want to process
    // the reads again and see if any discordant pairs support the split
    public void processReadGroup(@NotNull final ReadGroup rg)
    {
        // what we are looking for are reads and their mates

        TelomericBreakEnd tbe = findTelomericBreakEnd(r);

        if (tbe != null)
        {
            mPotentialBreakEnds.add(tbe);
            TE_LOGGER.trace("record: {} breakend: {}", r, tbe);
        }
    }

    // We want to check the read group against the break end see if this supports
    // any of the following scenarios:
    // for right C telomere
    public void checkAgainstBreakend(ReadGroup rg, TelomericBreakEnd breakEnd)
    {
        // firstly is any of the read mapped close to the break end?
        for (SAMRecord r : rg.Reads)
        {
            if (r.getReadUnmappedFlag())
            {
                continue;
            }
            // first check if the break end is contained within this read
            if (r.getAlignmentStart() >= breakEnd.getPosition() &&
                r.getAlignmentEnd() <= breakEnd.getPosition())
            {
                // this could be one of the reads that support this breakend in the first place
                // we need to filter it out
                continue;
            }

            if (Math.abs(r.getAlignmentStart() - breakEnd.getPosition()) < MAX_READ_DISTANCE ||
                Math.abs(r.getAlignmentEnd() - breakEnd.getPosition()) < MAX_READ_DISTANCE)
            {
                // we are kind of close to the break end

            }
        }

        // we are only interested in fragments that point
        // towards the telomere section from break end
        // fragments that point away from it won't be interesting
        switch (breakEnd.getType())
        {
            case RIGHT_G_TELOMERIC:
                // right side telomeric
                // we are interested in reads where the first in pair is before
                // the break end, and the 2nd in pair is after or somewhere else.

                break;
            case LEFT_C_TELOMERIC:
                break;
            case RIGHT_C_TELOMERIC:
                break;
            case LEFT_G_TELOMERIC:
                break;
        }
    }
}
