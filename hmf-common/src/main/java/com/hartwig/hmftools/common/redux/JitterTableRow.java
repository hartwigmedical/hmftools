package com.hartwig.hmftools.common.redux;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.common.bam.ConsensusType;

public class JitterTableRow
{
    private final int mRefNumUnits;
    private final String mRepeatUnit;
    private final ConsensusType mConsensusType;

    private int mTotalReadCount;
    private final Map<Integer,Integer> mJitterCounts;

    public JitterTableRow(final int refNumUnits, final String repeatUnit, final ConsensusType consensusType)
    {
        mRefNumUnits = refNumUnits;
        mRepeatUnit = repeatUnit;
        mConsensusType = consensusType;

        mTotalReadCount = 0;
        mJitterCounts = new HashMap<>();
    }

    public int refNumUnits() { return mRefNumUnits; }
    public String getRepeatUnit() { return mRepeatUnit; }
    public ConsensusType getConsensusType() { return mConsensusType; }
    public Map<Integer,Integer> jitterCounts() { return mJitterCounts; }

    public int totalReadCount() { return mTotalReadCount; }
    public void setTotalReadCount(int count) { mTotalReadCount = count; }

    public void addRead(int jitter)
    {
        addReads(jitter, 1);
    }

    public void addReads(int jitter, int numReads)
    {
        mTotalReadCount += numReads;
        mJitterCounts.merge(jitter, numReads, Integer::sum);
    }

    public int getJitterReadCount(int jitter)
    {
        return mJitterCounts.getOrDefault(jitter, 0);
    }

    public void setJitterReadCount(int jitter, int count)
    {
        mJitterCounts.put(jitter, count);
    }
}
