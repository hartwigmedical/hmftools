package com.hartwig.hmftools.common.basequal.jitter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

// works with one repeat
// find the number of reads of each number of repeat across the section
public class MicrosatelliteSiteAnalyser
{
    public static final Logger sLogger = LogManager.getLogger(MicrosatelliteSiteAnalyser.class);

    private final RefGenomeMicrosatellite mRefGenomeMicrosatellite;
    private final ConsensusMarker mConsensusMarker;
    private final List<MicrosatelliteRead> mMicrosatelliteReads;
    private final Set<String> mSeenConsensusTypes;

    public MicrosatelliteSiteAnalyser(final RefGenomeMicrosatellite refGenomeMicrosatellite,
            @Nullable final ConsensusMarker consensusMarker)
    {
        mRefGenomeMicrosatellite = refGenomeMicrosatellite;
        mConsensusMarker = consensusMarker;
        mMicrosatelliteReads = new ArrayList<>();
        mSeenConsensusTypes = Sets.newHashSet();
    }

    public RefGenomeMicrosatellite refGenomeMicrosatellite()
    {
        return mRefGenomeMicrosatellite;
    }

    public List<MicrosatelliteRead> getReadRepeatMatches()
    {
        return mMicrosatelliteReads;
    }

    public List<MicrosatelliteRead> getReadRepeatMatches(final ConsensusType consensusType)
    {
        return mMicrosatelliteReads.stream().filter(o -> o.consensusType() == consensusType).collect(Collectors.toList());
    }

    public List<MicrosatelliteRead> getPassingReadRepeatMatches()
    {
        return mMicrosatelliteReads.stream().filter(o -> !o.shouldDropRead).collect(Collectors.toList());
    }

    public List<MicrosatelliteRead> getPassingReadRepeatMatches(final ConsensusType consensusType)
    {
        return mMicrosatelliteReads.stream().filter(o -> !o.shouldDropRead && o.consensusType() == consensusType).collect(Collectors.toList());
    }

    public int numReadRejected()
    {
        return (int) mMicrosatelliteReads.stream().filter(o -> o.shouldDropRead).count();
    }

    public int numReadRejected(final ConsensusType consensusType)
    {
        return (int) mMicrosatelliteReads.stream().filter(o -> o.shouldDropRead && o.consensusType() == consensusType).count();
    }

    public Set<String> seenConsensusTypes() { return mSeenConsensusTypes; }

    public synchronized void addReadToStats(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag() || read.getDuplicateReadFlag())
        {
            return;
        }

        MicrosatelliteRead msRead = MicrosatelliteRead.from(mRefGenomeMicrosatellite, read, mConsensusMarker);
        if(msRead == null)
        {
            return;
        }

        mMicrosatelliteReads.add(msRead);
        if(!msRead.shouldDropRead)
        {
            mSeenConsensusTypes.add(msRead.consensusType().name());
        }
    }

    public int getCountWithRepeatUnits(int numRepeatUnits, ConsensusType consensusType)
    {
        return (int) getPassingReadRepeatMatches().stream().filter(o -> o.numRepeatUnits() == numRepeatUnits && o.consensusType() == consensusType).count();
    }

    public boolean shouldKeepSite(final double altCountFractionInit,
            final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff,
            final int minPassingSiteReads)
    {
        return !isRealVariant(altCountFractionInit, altCountFractionCutoffStep, rejectedReadFractionCutoff, minPassingSiteReads);
    }

    // have threshold for ALT site differ depending on INDEL length (e.g. 30% for INDEL=1, 25% for INDEL=2, ..., 10% for INDEL=5)
    public boolean isRealVariant(final double altCountFractionInit,
            final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff,
            final int minPassingSiteReads)
    {
        Validate.isTrue(altCountFractionCutoffStep <= 0.0);

        if(getPassingReadRepeatMatches().size() < minPassingSiteReads)
        {
            return true;
        }

        double fractionRejected = 1.0 - getPassingReadRepeatMatches().size() / (double) getReadRepeatMatches().size();

        if(Doubles.greaterOrEqual(fractionRejected, rejectedReadFractionCutoff))
        {
            return true;
        }

        Map<Integer, Integer> repeatReadCounts = new HashMap<>();

        for(MicrosatelliteRead microsatelliteRead : getPassingReadRepeatMatches())
        {
            int repeatDiff = mRefGenomeMicrosatellite.numRepeat - microsatelliteRead.numRepeatUnits();

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
