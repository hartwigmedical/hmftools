package com.hartwig.hmftools.common.basequal.jitter;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
    private final Set<String> mSeenConsensusTypes;

    private final Map<String, Map<Integer, Integer>> mPassingJitterCountsByConsensusType;
    private final Map<String, Integer> mReadRepeatMatchCountsByConsensusType;
    private final Map<String, Integer> mNumReadRejectedByConsensusType;
    private int mPassingReadRepeatMatchCount;

    private final Map<String, List<Integer>> mAllPassingRepeatLengthsByConsensusType;

    public MicrosatelliteSiteAnalyser(final RefGenomeMicrosatellite refGenomeMicrosatellite,
            @Nullable final ConsensusMarker consensusMarker, boolean storeAllPassingRepeatLengths)
    {
        mRefGenomeMicrosatellite = refGenomeMicrosatellite;
        mConsensusMarker = consensusMarker;
        mSeenConsensusTypes = Sets.newHashSet();

        mPassingJitterCountsByConsensusType = Maps.newHashMap();
        mReadRepeatMatchCountsByConsensusType = Maps.newHashMap();
        mNumReadRejectedByConsensusType = Maps.newHashMap();
        mPassingReadRepeatMatchCount = 0;

        mAllPassingRepeatLengthsByConsensusType = storeAllPassingRepeatLengths ? Maps.newHashMap() : null;
    }

    public RefGenomeMicrosatellite refGenomeMicrosatellite()
    {
        return mRefGenomeMicrosatellite;
    }

    public int readRepeatMatchCount(final ConsensusType consensusType)
    {
        return mReadRepeatMatchCountsByConsensusType.getOrDefault(consensusType.name(), 0);
    }

    public int readRepeatMatchCount()
    {
        return mReadRepeatMatchCountsByConsensusType.values().stream().mapToInt(x -> x).sum();
    }

    public int numReadRejected(final ConsensusType consensusType)
    {
        return mNumReadRejectedByConsensusType.getOrDefault(consensusType.name(), 0);
    }

    public Map<Integer, Integer> passingJitterCounts(final ConsensusType consensusType)
    {
        return mPassingJitterCountsByConsensusType.getOrDefault(consensusType.name(), Collections.emptyMap());
    }

    public int passingJitterCount(int jitter, final ConsensusType consensusType)
    {
        return passingJitterCounts(consensusType).getOrDefault(jitter, 0);
    }

    public Set<String> seenConsensusTypes() { return mSeenConsensusTypes; }

    public List<Integer> allPassingRepeats(final ConsensusType consensusType)
    {
        return mAllPassingRepeatLengthsByConsensusType.getOrDefault(consensusType.name(), Collections.emptyList());
    }

    public synchronized void addReadToStats(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag() || read.getDuplicateReadFlag())
            return;

        MicrosatelliteRead msRead = MicrosatelliteRead.from(mRefGenomeMicrosatellite, read, mConsensusMarker);
        if(msRead == null)
            return;

        ConsensusType consensusType = msRead.consensusType();
        mReadRepeatMatchCountsByConsensusType.merge(consensusType.name(), 1, Integer::sum);

        if(msRead.shouldDropRead)
        {
            mNumReadRejectedByConsensusType.merge(consensusType.name(), 1, Integer::sum);
            return;
        }

        mSeenConsensusTypes.add(msRead.consensusType().name());
        mPassingReadRepeatMatchCount++;
        mPassingJitterCountsByConsensusType.computeIfAbsent(consensusType.name(), key -> Maps.newHashMap());
        mPassingJitterCountsByConsensusType.get(consensusType.name()).merge(msRead.jitter(), 1, Integer::sum);

        if(mAllPassingRepeatLengthsByConsensusType != null)
        {
            mAllPassingRepeatLengthsByConsensusType.computeIfAbsent(consensusType.name(), key -> Lists.newArrayList());
            mAllPassingRepeatLengthsByConsensusType.get(consensusType.name()).add(msRead.readRepeatLength());
        }
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

        if(mPassingReadRepeatMatchCount < minPassingSiteReads)
        {
            return true;
        }

        double fractionRejected = 1.0 - mPassingReadRepeatMatchCount / (double) readRepeatMatchCount();

        if(Doubles.greaterOrEqual(fractionRejected, rejectedReadFractionCutoff))
        {
            return true;
        }

        Map<Integer, Integer> passingJitterCounts = Maps.newHashMap();
        for(Map<Integer, Integer> jitterCounts : mPassingJitterCountsByConsensusType.values())
        {
            for(Map.Entry<Integer, Integer> entry : jitterCounts.entrySet())
            {
                int jitter = entry.getKey();
                if(jitter == 0)
                    continue;

                int numReads = entry.getValue();
                passingJitterCounts.merge(jitter, numReads, Integer::sum);
            }
        }

        for(Map.Entry<Integer, Integer> entry : passingJitterCounts.entrySet())
        {
            int repeatDiff = entry.getKey();
            int readCount = entry.getValue();

            double fractionCutoff = Math.max(altCountFractionInit + (Math.abs(repeatDiff) - 1) * altCountFractionCutoffStep, 0.1);
            double countCutoff = fractionCutoff * mPassingReadRepeatMatchCount;
            if(Doubles.greaterThan(readCount, countCutoff))
            {
                return true;
            }
        }

        return false;
    }
}
