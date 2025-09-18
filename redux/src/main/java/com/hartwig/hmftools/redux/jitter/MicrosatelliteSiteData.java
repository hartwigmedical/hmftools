package com.hartwig.hmftools.redux.jitter;

import java.util.Collections;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.utils.Doubles;

import htsjdk.samtools.SAMRecord;

public class MicrosatelliteSiteAnalyser
{
    private final MicrosatelliteSite mMicrosatelliteSite;
    private final ConsensusMarker mConsensusMarker;

    private final EnumMap<ConsensusType, Map<Integer,Integer>> mPassingJitterCountsByConsensusType;
    private final EnumMap<ConsensusType, Integer> mReadRepeatMatchCountsByConsensusType;
    private int mNumReadRejected;
    private int mPassingReadRepeatMatchCount;

    private final EnumMap<ConsensusType, SortedMap<Integer, Integer>> mPassingRepeatLengthCountsByConsensusType;

    public MicrosatelliteSiteAnalyser(
            final MicrosatelliteSite microsatelliteSite, final ConsensusMarker consensusMarker,
            boolean storeAllPassingRepeatLengths)
    {
        mMicrosatelliteSite = microsatelliteSite;
        mConsensusMarker = consensusMarker;

        mPassingJitterCountsByConsensusType = Maps.newEnumMap(ConsensusType.class);
        mReadRepeatMatchCountsByConsensusType = Maps.newEnumMap(ConsensusType.class);
        mNumReadRejected = 0;
        mPassingReadRepeatMatchCount = 0;

        mPassingRepeatLengthCountsByConsensusType = storeAllPassingRepeatLengths ? Maps.newEnumMap(ConsensusType.class) : null;
    }

    public MicrosatelliteSite refGenomeMicrosatellite()
    {
        return mMicrosatelliteSite;
    }

    public int readRepeatMatchCount(final ConsensusType consensusType)
    {
        return mReadRepeatMatchCountsByConsensusType.getOrDefault(consensusType, 0);
    }

    public int readRepeatMatchCount()
    {
        return mReadRepeatMatchCountsByConsensusType.values().stream().mapToInt(x -> x).sum();
    }

    public int numReadRejected() { return mNumReadRejected; }

    public Map<Integer,Integer> passingJitterCounts(final ConsensusType consensusType)
    {
        return mPassingJitterCountsByConsensusType.getOrDefault(consensusType, Collections.emptyMap());
    }

    public int passingJitterCount(int jitter, final ConsensusType consensusType)
    {
        return passingJitterCounts(consensusType).getOrDefault(jitter, 0);
    }

    public List<Integer> allPassingRepeats(final ConsensusType consensusType)
    {
        if(!mPassingRepeatLengthCountsByConsensusType.containsKey(consensusType))
            return Collections.emptyList();

        List<Integer> passingRepeats = Lists.newArrayList();

        for(Map.Entry<Integer, Integer> repeatLengthAndCount : mPassingRepeatLengthCountsByConsensusType.get(consensusType).entrySet())
        {
            int repeatLength = repeatLengthAndCount.getKey();
            int repeatCount = repeatLengthAndCount.getValue();
            for(int i = 0; i < repeatCount; i++)
            {
                passingRepeats.add(repeatLength);
            }
        }

        return passingRepeats;
    }

    public synchronized void addReadToStats(final SAMRecord read)
    {
        MicrosatelliteRead msRead = MicrosatelliteRead.from(mMicrosatelliteSite, read, mConsensusMarker);

        ConsensusType consensusType = msRead.consensusType();
        mReadRepeatMatchCountsByConsensusType.merge(consensusType, 1, Integer::sum);

        if(!msRead.isValidRead())
        {
            ++mNumReadRejected;
            return;
        }

        mPassingReadRepeatMatchCount++;
        mPassingJitterCountsByConsensusType.computeIfAbsent(consensusType, key -> Maps.newHashMap());
        mPassingJitterCountsByConsensusType.get(consensusType).merge(msRead.jitter(), 1, Integer::sum);

        if(mPassingRepeatLengthCountsByConsensusType != null)
        {
            int readRepeatLength = msRead.readRepeatLength();
            mPassingRepeatLengthCountsByConsensusType.computeIfAbsent(consensusType, key -> Maps.newTreeMap());
            mPassingRepeatLengthCountsByConsensusType.get(consensusType).merge(readRepeatLength, 1, Integer::sum);
        }
    }

    public boolean shouldKeepSite(
            final double altCountFractionInit, final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff, final int minPassingSiteReads)
    {
        return !isRealVariant(altCountFractionInit, altCountFractionCutoffStep, rejectedReadFractionCutoff, minPassingSiteReads);
    }

    public boolean isRealVariant(
            final double altCountFractionInit, final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff, final int minPassingSiteReads)
    {
        // have threshold for ALT site differ depending on INDEL length (e.g. 30% for INDEL=1, 25% for INDEL=2, ..., 10% for INDEL=5)
        if(mPassingReadRepeatMatchCount < minPassingSiteReads)
            return true;

        double fractionRejected = 1.0 - mPassingReadRepeatMatchCount / (double) readRepeatMatchCount();

        if(Doubles.greaterOrEqual(fractionRejected, rejectedReadFractionCutoff))
            return true;

        Map<Integer,Integer> passingJitterCounts = Maps.newHashMap();

        for(Map<Integer,Integer> jitterCounts : mPassingJitterCountsByConsensusType.values())
        {
            for(Map.Entry<Integer,Integer> entry : jitterCounts.entrySet())
            {
                int jitter = entry.getKey();
                if(jitter == 0)
                    continue;

                int numReads = entry.getValue();
                passingJitterCounts.merge(jitter, numReads, Integer::sum);
            }
        }

        for(Map.Entry<Integer,Integer> entry : passingJitterCounts.entrySet())
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
