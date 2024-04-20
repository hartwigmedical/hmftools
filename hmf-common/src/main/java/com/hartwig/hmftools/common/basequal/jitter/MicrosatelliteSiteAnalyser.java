package com.hartwig.hmftools.common.basequal.jitter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;

// works with one repeat
// find the number of reads of each number of repeat across the section
public class MicrosatelliteSiteAnalyser
{
    public static final Logger sLogger = LogManager.getLogger(MicrosatelliteSiteAnalyser.class);

    final RefGenomeMicrosatellite refGenomeMicrosatellite;

    private final List<MicrosatelliteRead> mMicrosatelliteReads = new ArrayList<>();

    public List<MicrosatelliteRead> getReadRepeatMatches() { return mMicrosatelliteReads; }

    public List<MicrosatelliteRead> getPassingReadRepeatMatches()
    {
        return mMicrosatelliteReads.stream().filter(o -> !o.shouldDropRead).collect(Collectors.toList());
    }

    public int numReadRejected()
    {
        return (int) mMicrosatelliteReads.stream().filter(o -> o.shouldDropRead).count();
    }

    public MicrosatelliteSiteAnalyser(final RefGenomeMicrosatellite refGenomeMicrosatellite)
    {
        this.refGenomeMicrosatellite = refGenomeMicrosatellite;
    }

    public synchronized void addReadToStats(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag() || read.getDuplicateReadFlag())
            return;

        mMicrosatelliteReads.add(MicrosatelliteRead.from(refGenomeMicrosatellite, read));
    }

    public int getCountWithRepeatUnits(int numRepeatUnits)
    {
        return (int)getPassingReadRepeatMatches().stream().filter(o -> o.numRepeatUnits() == numRepeatUnits).count();
    }

    public boolean shouldKeepSite(final double altCountFractionInit,
            final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff)
    {
        return !isRealVariant(altCountFractionInit, altCountFractionCutoffStep, rejectedReadFractionCutoff);
    }

    // have threshold for ALT site differ depending on INDEL length (e.g. 30% for INDEL=1, 25% for INDEL=2, ..., 10% for INDEL=5)
    public boolean isRealVariant(final double altCountFractionInit,
            final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff)
    {
        Validate.isTrue(altCountFractionCutoffStep <= 0.0);

        double fractionRejected = 1.0 - getPassingReadRepeatMatches().size() / (double)getReadRepeatMatches().size();

        if(Doubles.greaterOrEqual(fractionRejected, rejectedReadFractionCutoff))
        {
            return true;
        }

        Map<Integer, Integer> repeatReadCounts = new HashMap<>();

        for(MicrosatelliteRead microsatelliteRead : getPassingReadRepeatMatches())
        {
            int repeatDiff = refGenomeMicrosatellite.numRepeat - microsatelliteRead.numRepeatUnits();

            if(repeatDiff != 0)
            {
                repeatReadCounts.merge(repeatDiff, 1, Integer::sum);
            }
        }

        for(Map.Entry<Integer, Integer> entry : repeatReadCounts.entrySet())
        {
            int repeatDiff = entry.getKey();
            int readCount = entry.getValue();

            double fractionCutoff = Math.max(altCountFractionInit + (Math.abs(repeatDiff) - 1) * altCountFractionCutoffStep, 0.1);
            double countCutoff = fractionCutoff * getPassingReadRepeatMatches().size();
            if(Doubles.greaterThan(readCount, countCutoff))
            {
                return true;
            }
        }

        return false;
    }
}
