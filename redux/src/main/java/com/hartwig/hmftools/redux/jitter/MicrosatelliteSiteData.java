package com.hartwig.hmftools.redux.jitter;

import static com.hartwig.hmftools.redux.ReduxConfig.isSbx;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.ALT_COUNT_FRACTION_INIT;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.ALT_COUNT_FRACTION_STEP;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.MAX_REJECTED_READ_FRACTION;
import static com.hartwig.hmftools.redux.jitter.JitterConstants.MIN_PASSING_SITE_READS;

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

public class MicrosatelliteSiteData
{
    private final MicrosatelliteSite mMicrosatelliteSite;
    private final ConsensusMarker mConsensusMarker;

    private int mTotalReadCount;
    private int mPassingReadCount;

    private final EnumMap<ConsensusType,Map<Integer,Integer>> mJitterCountsByConsensusType;
    private final EnumMap<ConsensusType,Integer> mCountsByConsensusType;
    private final EnumMap<ConsensusType,Integer> mRejectedCountsByConsensusType;

    private final EnumMap<ConsensusType, SortedMap<Integer, Integer>> mRepeatLengthCountsByConsensusType;

    public MicrosatelliteSiteData(
            final MicrosatelliteSite microsatelliteSite, final ConsensusMarker consensusMarker,
            boolean storeAllPassingRepeatLengths)
    {
        mMicrosatelliteSite = microsatelliteSite;
        mConsensusMarker = consensusMarker;

        mTotalReadCount = 0;
        mPassingReadCount = 0;

        mJitterCountsByConsensusType = Maps.newEnumMap(ConsensusType.class);
        mCountsByConsensusType = Maps.newEnumMap(ConsensusType.class);
        mRejectedCountsByConsensusType = Maps.newEnumMap(ConsensusType.class);

        mRepeatLengthCountsByConsensusType = storeAllPassingRepeatLengths ? Maps.newEnumMap(ConsensusType.class) : null;
    }

    public MicrosatelliteSite refGenomeMicrosatellite() { return mMicrosatelliteSite; }

    public int readCountByConsensus(final ConsensusType consensusType)
    {
        return mCountsByConsensusType.getOrDefault(consensusType, 0);
    }

    public int totalReadCount()
    {
        return mTotalReadCount;
    }

    public int readsRejectedByConsensus(final ConsensusType consensusType)
    {
        return mRejectedCountsByConsensusType.getOrDefault(consensusType, 0);
    }

    public Map<Integer,Integer> passingJitterCounts(final ConsensusType consensusType)
    {
        return mJitterCountsByConsensusType.getOrDefault(consensusType, Collections.emptyMap());
    }

    public int passingJitterCount(int jitter, final ConsensusType consensusType)
    {
        return passingJitterCounts(consensusType).getOrDefault(jitter, 0);
    }

    public List<Integer> allPassingRepeats(final ConsensusType consensusType)
    {
        if(!mRepeatLengthCountsByConsensusType.containsKey(consensusType))
            return Collections.emptyList();

        List<Integer> passingRepeats = Lists.newArrayList();

        for(Map.Entry<Integer, Integer> repeatLengthAndCount : mRepeatLengthCountsByConsensusType.get(consensusType).entrySet())
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

        if(!msRead.isValidRead())
        {
            Integer rejectedCount = mRejectedCountsByConsensusType.get(consensusType);
            mRejectedCountsByConsensusType.put(consensusType, rejectedCount != null ? rejectedCount + 1 : 1);
            return;
        }

        mPassingReadCount++;

        if(isSbx())
        {
            // count towards passing count but nothing else if formed from simplex reads
            if(consensusType == ConsensusType.SINGLE || msRead.hasMediumQualBases())
                return;
        }

        Integer consensusCount = mCountsByConsensusType.get(consensusType);
        mCountsByConsensusType.put(consensusType, consensusCount != null ? consensusCount + 1 : 1);

        Map<Integer,Integer> consensusLengthCounts = mJitterCountsByConsensusType.get(consensusType);

        if(consensusLengthCounts == null)
        {
            consensusLengthCounts = Maps.newHashMap();
            mJitterCountsByConsensusType.put(consensusType, consensusLengthCounts);
        }

        Integer lengthCount = consensusLengthCounts.get(msRead.jitter());
        consensusLengthCounts.put(msRead.jitter(), lengthCount != null ? lengthCount + 1 : 1);

        if(mRepeatLengthCountsByConsensusType != null)
        {
            int readRepeatLength = msRead.readRepeatLength();
            mRepeatLengthCountsByConsensusType.computeIfAbsent(consensusType, key -> Maps.newTreeMap());
            mRepeatLengthCountsByConsensusType.get(consensusType).merge(readRepeatLength, 1, Integer::sum);
        }
    }

    public boolean shouldKeepSite()
    {
        // return !isRealVariant(altCountFractionInit, altCountFractionCutoffStep, rejectedReadFractionCutoff, minPassingSiteReads);
        return !isRealVariant(ALT_COUNT_FRACTION_INIT, ALT_COUNT_FRACTION_STEP, MAX_REJECTED_READ_FRACTION, MIN_PASSING_SITE_READS);
    }

    public boolean isRealVariant(
            final double altCountFractionInit, final double altCountFractionCutoffStep,
            final double rejectedReadFractionCutoff, final int minPassingSiteReads)
    {
        // have threshold for ALT site differ depending on INDEL length (e.g. 30% for INDEL=1, 25% for INDEL=2, ..., 10% for INDEL=5)
        if(mPassingReadCount < minPassingSiteReads)
            return true;

        double fractionRejected = 1.0 - mPassingReadCount / (double)mTotalReadCount;

        if(Doubles.greaterOrEqual(fractionRejected, rejectedReadFractionCutoff))
            return true;

        Map<Integer,Integer> passingJitterCounts = Maps.newHashMap();

        for(Map<Integer,Integer> jitterCounts : mJitterCountsByConsensusType.values())
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
            double countCutoff = fractionCutoff * mPassingReadCount;
            if(Doubles.greaterThan(readCount, countCutoff))
            {
                return true;
            }
        }

        return false;
    }
}
