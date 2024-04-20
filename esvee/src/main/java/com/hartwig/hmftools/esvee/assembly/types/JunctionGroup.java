package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConstants.BAM_READ_JUNCTION_BUFFER;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class JunctionGroup extends BaseRegion
{
    private final List<Junction> mJunctions;
    private final List<Read> mCandidateReads;

    private int mIndex; // in the full list for the chromosome

    private final List<JunctionAssembly> mJunctionAssemblies;
    private final List<DiscordantGroup> mDiscordantGroups;

    public JunctionGroup(final Junction junction)
    {
        super(junction.Position, junction.Position);
        mJunctions = Lists.newArrayList(junction);
        mCandidateReads = Lists.newArrayList();
        mJunctionAssemblies = Lists.newArrayList();
        mDiscordantGroups = Lists.newArrayList();
        mIndex = -1;
    }

    public List<Junction> junctions() { return mJunctions; }
    public int count() { return mJunctions.size(); }

    public String chromosome() { return mJunctions.get(0).Chromosome; }

    public int minPosition() { return start(); }
    public int maxPosition() { return end(); }
    public int range() { return super.baseLength(); }

    public int readRangeStart() { return start() - BAM_READ_JUNCTION_BUFFER; }
    public int readRangeEnd() { return end() + BAM_READ_JUNCTION_BUFFER; }

    public void setIndex(int index) { mIndex = index; }
    public int index() { return mIndex; }

    public void addJunction(final Junction junction)
    {
        mJunctions.add(junction);
        setEnd(max(end(), junction.Position));
        setStart(min(start(), junction.Position));
    }

    public List<Read> candidateReads() { return mCandidateReads; }
    public void addCandidateRead(final Read read) { mCandidateReads.add(read); }
    public void clearCandidateReads() { mCandidateReads.clear(); }
    public int candidateReadCount() { return mCandidateReads.size(); }

    public void addJunctionAssemblies(final List<JunctionAssembly> assemblies) { mJunctionAssemblies.addAll(assemblies); }
    public List<JunctionAssembly> junctionAssemblies() { return mJunctionAssemblies; }

    public void addDiscordantGroups(final List<DiscordantGroup> groups) { mDiscordantGroups.addAll(groups); }
    public List<DiscordantGroup> discordantGroups() { return mDiscordantGroups; }

    public boolean overlapsRemoteRegion(final RemoteRegion region)
    {
        return positionsOverlap(readRangeStart(), readRangeEnd(), region.start(), region.end());
    }

    public String toString()
    {
        return format("%s:%d-%d range(%d) count(%d) assemblies(%d) reads(%d)",
            chromosome(), start(), end(), range(), mJunctions.size(), mJunctionAssemblies.size(), mCandidateReads.size());
    }

    public static <E extends BaseRegion> int binarySearch(int readStart, final List<E> regions)
    {
        // Returns index of the last region in regions with start pos less than or equal to readStart
        // If all regions have start pos larger than readStart then return zero. It is assumed that regions are sorted
        int binarySearchIndex = Collections.binarySearch(regions, new BaseRegion(readStart, readStart));

        if(binarySearchIndex >= 0)
            return binarySearchIndex; // found with exact match for start pos

        // get insertion point
        int insertionIndex = -(binarySearchIndex + 1);

        if(insertionIndex == 0)
            return 0;

        return insertionIndex - 1;
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
