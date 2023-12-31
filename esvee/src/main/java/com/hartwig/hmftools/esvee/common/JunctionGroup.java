package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

import org.jetbrains.annotations.NotNull;

public class JunctionGroup implements Comparable<JunctionGroup>
{
    private final List<Junction> mJunctions;

    private int mMinPosition;
    private int mMaxPosition;

    private final List<Read> mCandidateReads;

    public JunctionGroup(final Junction junction)
    {
        mJunctions = Lists.newArrayList(junction);
        mMaxPosition = junction.Position;
        mMinPosition = junction.Position;
        mCandidateReads = Lists.newArrayList();
    }

    public List<Junction> junctions() { return mJunctions; }
    public int count() { return mJunctions.size(); }

    public String chromosome() { return mJunctions.get(0).Chromosome; }

    public int minPosition() { return mMinPosition; }
    public int maxPosition() { return mMaxPosition; }
    public int range() { return mMaxPosition - mMinPosition + 1; }

    public void addJunction(final Junction junction)
    {
        mJunctions.add(junction);
        mMaxPosition = max(mMaxPosition, junction.Position);
        mMinPosition = min(mMinPosition, junction.Position);
    }

    public List<Read> candidateReads() { return mCandidateReads; }
    public void addCandidateRead(final Read read) { mCandidateReads.add(read); }
    public void clearCandidateReads() { mCandidateReads.clear(); }
    public int candidateReadCount() { return mCandidateReads.size(); }

    public String toString()
    {
        return format("%s:%d-%d range(%d) count(%d) reads(%d)",
            chromosome(), mMinPosition, mMaxPosition, range(), mJunctions.size(), mCandidateReads.size());
    }

    @Override
    public int compareTo(@NotNull final JunctionGroup other)
    {
        int firstRange = range();
        int secondRange = other.range();

        if(firstRange == secondRange)
            return 0;

        return firstRange > secondRange ? -1 : 1;
    }

    public static List<JunctionGroup> buildJunctionGroups(final List<Junction> junctions, final int maxDistance)
    {
        List<JunctionGroup> groups = Lists.newArrayList();

        JunctionGroup currentGroup = new JunctionGroup(junctions.get(0));
        groups.add(currentGroup);

        for(int i = 1; i < junctions.size(); ++i)
        {
            Junction junction = junctions.get(i);

            if(junction.Position - currentGroup.maxPosition() <= maxDistance)
            {
                currentGroup.addJunction(junction);
            }
            else
            {
                currentGroup = new JunctionGroup(junction);
                groups.add(currentGroup);
            }
        }

        return groups;
    }
}
