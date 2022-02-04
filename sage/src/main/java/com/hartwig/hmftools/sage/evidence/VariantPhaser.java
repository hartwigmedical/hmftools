package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;

public class VariantPhaser
{
    private final PhaseSetCounter mPhaseSetCounter;

    private boolean mEnabled;
    private boolean mLogData;
    private ChrBaseRegion mRegion; // only for logging

    private final List<PhasedReadCounters> mPhasedReadCounters; // order by position of the lowest variant
    private int mCurrentIndex;

    public static final int MIN_READ_COUNT = 1;

    public VariantPhaser(final PhaseSetCounter phaseSetCounter)
    {
        mPhaseSetCounter = phaseSetCounter;
        mEnabled = false;
        mLogData = false;

        mPhasedReadCounters = Lists.newArrayList();
        mCurrentIndex = 0;
    }

    public boolean enabled() { return mEnabled; }

    public void initialise(final ChrBaseRegion region, boolean enabled, boolean logData)
    {
        mRegion = region;
        mEnabled = enabled;
        mLogData = logData;
        mPhasedReadCounters.clear();
        mCurrentIndex = 0;
    }

    public void registeredPhasedVariants(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        if(!mEnabled)
            return;

        if(posCounters.isEmpty() || posCounters.size() + negCounters.size() < 2)
            return;

        // walk backwards then forwards from the current location looking for a match
        int minVarPos = posCounters.get(0).position();
        int maxVarPos = posCounters.get(posCounters.size() - 1).position();

        if(!negCounters.isEmpty())
        {
            minVarPos = min(minVarPos, posCounters.get(0).position());
            maxVarPos = max(maxVarPos, negCounters.get(negCounters.size() - 1).position());
        }

        if(mPhasedReadCounters.isEmpty())
        {
            mPhasedReadCounters.add(new PhasedReadCounters(mPhasedReadCounters.size(), minVarPos, maxVarPos, posCounters, negCounters));
            return;
        }

        int index = mCurrentIndex;

        PhasedReadCounters currentPhasedRc = mPhasedReadCounters.get(index);

        boolean searchBackwards = minVarPos < currentPhasedRc.MinVariantPos;

        boolean matchFound = false;

        if(searchBackwards)
        {
            while(index >= 0)
            {
                PhasedReadCounters phasedRc = mPhasedReadCounters.get(index);

                if(phasedRc.MinVariantPos < minVarPos)
                {
                    ++index;
                    break;
                }

                if(phasedRc.MinVariantPos == minVarPos)
                {
                    // test for an exact match
                    if(phasedRc.isExactMatch(minVarPos, maxVarPos, posCounters, negCounters))
                    {
                        phasedRc.ReadCount++;
                        matchFound = true;
                        break;
                    }
                }

                index--;

                if(index < 0)
                {
                    index = 0;
                    break;
                }
            }
        }
        else
        {
            while(index < mPhasedReadCounters.size())
            {
                PhasedReadCounters phasedRc = mPhasedReadCounters.get(index);

                if(minVarPos < phasedRc.MinVariantPos)
                    break;

                if(phasedRc.MinVariantPos == minVarPos)
                {
                    // test for an exact match
                    if(phasedRc.isExactMatch(minVarPos, maxVarPos, posCounters, negCounters))
                    {
                        phasedRc.ReadCount++;
                        matchFound = true;
                        break;
                    }
                }

                index++;
            }
        }

        mCurrentIndex = index;

        if(matchFound)
            return;

        if(index > 0 && index < mPhasedReadCounters.size() - 1)
        {
            if(minVarPos < mPhasedReadCounters.get(index - 1).MinVariantPos || minVarPos > mPhasedReadCounters.get(index + 1).MinVariantPos)
            {
                SG_LOGGER.error("invalid insertion");
            }
        }

        mPhasedReadCounters.add(index, new PhasedReadCounters(mPhasedReadCounters.size(), minVarPos, maxVarPos, posCounters, negCounters));
    }

    private void mergeGroups()
    {
        /* merge groups by the following rules:
            - if a group has only matching positive read-counters (ie not a subset of any other group), then merge them all
            - if a group is a subset of 1 or more other groups, then allocate its reads pro-rata and remove it
            - finally merge any pair which share a subset (without either being a subset of each other)

        then:
            - filter for LPS support > X
            - remove any uninformative LPS (ie LPS has 1+ variant and is the only LPS that includes that variant).
        */

        Set<PhasedReadCounters> removedGroups = Sets.newHashSet();

        /*
        List<PhasedReadCounters> phasedReadCounters = mPhasedReadCountersMap.values().stream().collect(Collectors.toList());
        Collections.sort(phasedReadCounters, new PhasedRcReadCountComparator());
        Collections.sort(phasedReadCounters, new PhasedRcPosCountComparator());
        Collections.sort(phasedReadCounters, new PhasedRcNegCountComparator());
        */

        for(PhasedReadCounters phasedRc : mPhasedReadCounters)
        {
            if(removedGroups.contains(phasedRc))
                continue;

            // check for super-set groups
            List<PhasedReadCounters> superGroups = mPhasedReadCounters.stream()
                    .filter(x -> x != phasedRc)
                    .filter(x -> !removedGroups.contains(x))
                    .filter(x -> phasedRc.positionsOverlap(x))
                    .filter(x -> phasedRc.isSubsetOf(x))
                    .collect(Collectors.toList());

            if(superGroups.isEmpty())
                continue;

            // where the super groups differ, allocate by pro-rata and continue
            if(superGroups.stream().anyMatch(x -> x.PositiveReadCounters.size() > phasedRc.PositiveReadCounters.size()))
            {
                double totalReads = superGroups.stream().mapToInt(x -> x.ReadCount).sum();

                for(PhasedReadCounters superGroup : superGroups)
                {
                    double allocFraction = superGroup.ReadCount / totalReads;
                    superGroup.AllocatedReadCount += allocFraction * phasedRc.ReadCount;
                    superGroup.AllocatedReadCount += allocFraction * phasedRc.AllocatedReadCount;
                }

                // remove subset group
                removedGroups.add(phasedRc);
            }
            else
            {
                // otherwise collapse the groups into this one
                for(PhasedReadCounters superGroup : superGroups)
                {
                    phasedRc.ReadCount += superGroup.ReadCount;
                    phasedRc.AllocatedReadCount += superGroup.AllocatedReadCount;
                }

                removedGroups.addAll(superGroups);
            }
        }

        removedGroups.forEach(x -> mPhasedReadCounters.remove(x));

        int i = 0;
        while(i < mPhasedReadCounters.size())
        {
            PhasedReadCounters phasedRc = mPhasedReadCounters.get(i);

            // remove any group or those with a single variant not present as a +ve in another group
            if(phasedRc.PositiveReadCounters.size() == 1)
            {
                if(mPhasedReadCounters.stream()
                        .filter(x -> x != phasedRc)
                        .filter(x -> phasedRc.positionsOverlap(x))
                        .noneMatch(x -> x.PositiveReadCounters.contains(phasedRc)))
                {
                    mPhasedReadCounters.remove(i);
                    continue;
                }
            }

            int j = i + 1;
            while(j < mPhasedReadCounters.size())
            {
                PhasedReadCounters otherPhasedRc = mPhasedReadCounters.get(j);

                if(!phasedRc.positionsOverlap(otherPhasedRc))
                    break;

                if(phasedRc.haveCommonSubset(otherPhasedRc))
                {
                    phasedRc.merge(otherPhasedRc);
                    mPhasedReadCounters.remove(j);
                }
                else
                {
                    ++j;
                }
            }

            ++i;
        }
    }

    public void assignLocalPhaseSets()
    {
        // assign local phase set IDs to all phased variants
        if(!mEnabled)
            return;

        int startCount = mPhasedReadCounters.size();

        mergeGroups();

        if(mLogData)
            logPhasedReadCounters();

        int assignedLps = 0;

        for(PhasedReadCounters phasedRc : mPhasedReadCounters)
        {
            if(phasedRc.ReadCount + phasedRc.AllocatedReadCount < MIN_READ_COUNT)
                continue;

            int nextLps = mPhaseSetCounter.getNext();
            phasedRc.PositiveReadCounters.forEach(x -> x.addLocalPhaseSet(nextLps, phasedRc.ReadCount, phasedRc.AllocatedReadCount));
            ++assignedLps;
        }

        SG_LOGGER.trace("region({}) phasing groups start({}) postMerge({}) assigned({})",
                mRegion, startCount, mPhasedReadCounters.size(), assignedLps);

        mPhasedReadCounters.clear();
    }

    private class PhasedReadCounters
    {
        public final int Id;
        public final int MinVariantPos;
        public final int MaxVariantPos;
        public final List<ReadContextCounter> PositiveReadCounters; // supported by the reads
        public final List<ReadContextCounter> NegativeReadCounters; // not supported by the reads

        public int ReadCount; // from uniquely supporting reads
        public double AllocatedReadCount; // allocated from subset groups

        public PhasedReadCounters(
                final int id, final int minVariantPos, final int maxVariantPos,
                final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
        {
            Id = id;
            MinVariantPos = minVariantPos;
            MaxVariantPos = maxVariantPos;
            PositiveReadCounters = posCounters;
            NegativeReadCounters = negCounters;
            ReadCount = 1;
            AllocatedReadCount = 0;
        }

        public boolean positionsOverlap(final PhasedReadCounters other)
        {
            return BaseRegion.positionsOverlap(MinVariantPos, MaxVariantPos, other.MinVariantPos, other.MaxVariantPos);
        }

        public boolean isExactMatch(
                final int minVariantPos, final int maxVariantPos,
                final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
        {
            if(minVariantPos != MinVariantPos || maxVariantPos != MaxVariantPos)
                return false;

            if(PositiveReadCounters.size() != posCounters.size() || NegativeReadCounters.size() != negCounters.size())
                return false;

            if(PositiveReadCounters.stream().anyMatch(x -> !posCounters.contains(x)))
                return false;

            if(NegativeReadCounters.stream().anyMatch(x -> !negCounters.contains(x)))
                return false;

            return true;
        }

        public void merge(final PhasedReadCounters other)
        {
            ReadCount += other.ReadCount;
            AllocatedReadCount += other.AllocatedReadCount;

            other.PositiveReadCounters.stream().filter(x -> !PositiveReadCounters.contains(x)).forEach(x -> PositiveReadCounters.add(x));
            other.NegativeReadCounters.stream().filter(x -> !NegativeReadCounters.contains(x)).forEach(x -> NegativeReadCounters.add(x));

            // leave coords as-is since they're not used during or after merging
        }

        private boolean isSubsetOf(final PhasedReadCounters other)
        {
            // returns true if this group is a subset of 'other' is a su
            if(other.PositiveReadCounters.size() < PositiveReadCounters.size())
                return false;

            if(!PositiveReadCounters.stream().allMatch(x -> other.PositiveReadCounters.contains(x)))
                return false;

            // cannot have contradictory negatives
            if(hasAnyOverlap(other.PositiveReadCounters, NegativeReadCounters) || hasAnyOverlap(PositiveReadCounters, other.NegativeReadCounters))
                return false;

            return true;
        }

        private boolean haveCommonSubset(final PhasedReadCounters other)
        {
            // returns true if the groups contain common subsets but not all
            int countInOther = (int)PositiveReadCounters.stream().filter(x -> other.PositiveReadCounters.contains(x)).count();

            if(countInOther == 0 || countInOther == PositiveReadCounters.size() || countInOther == other.PositiveReadCounters.size())
                return false;

            // cannot have contradictory negatives
            if(hasAnyOverlap(other.PositiveReadCounters, NegativeReadCounters) || hasAnyOverlap(PositiveReadCounters, other.NegativeReadCounters))
                return false;

            return true;
        }

        private boolean hasAnyOverlap(final List<ReadContextCounter> counters1, final List<ReadContextCounter> counters2)
        {
            return counters1.stream().anyMatch(x -> counters2.contains(x)) || counters2.stream().anyMatch(x -> counters1.contains(x));
        }

        public String toString()
        {
            return String.format("%d: range(%d - %d) pos(%d) neg(%d) rc(%d) alloc(%.1f)",
                Id, MinVariantPos, MaxVariantPos, PositiveReadCounters.size(), NegativeReadCounters.size(), ReadCount, AllocatedReadCount);
        }
    }

    public static class ReadCounterIdComparator implements Comparator<ReadContextCounter>
    {
        public int compare(final ReadContextCounter first, final ReadContextCounter second)
        {
            return first.id() - second.id();
        }
    }

    public static class PhasedRcReadCountComparator implements Comparator<PhasedReadCounters>
    {
        public int compare(final PhasedReadCounters first, final PhasedReadCounters second)
        {
            return second.ReadCount - first.ReadCount;
        }
    }

    public static class PhasedRcPosCountComparator implements Comparator<PhasedReadCounters>
    {
        public int compare(final PhasedReadCounters first, final PhasedReadCounters second)
        {
            return second.PositiveReadCounters.size() - first.PositiveReadCounters.size();
        }
    }

    public static class PhasedRcNegCountComparator implements Comparator<PhasedReadCounters>
    {
        public int compare(final PhasedReadCounters first, final PhasedReadCounters second)
        {
            return second.NegativeReadCounters.size() - first.NegativeReadCounters.size();
        }
    }

    private void logPhasedReadCounters()
    {
        Map<ReadContextCounter,Integer> rcIds = Maps.newHashMap();

        for(PhasedReadCounters phasedReadCounters : mPhasedReadCounters)
        {
            StringJoiner posVars = new StringJoiner(";");
            StringJoiner posIds = new StringJoiner(";");
            StringJoiner negVars = new StringJoiner(";");
            StringJoiner negIds = new StringJoiner(";");

            for(ReadContextCounter rc : phasedReadCounters.PositiveReadCounters)
            {
                posVars.add(getReadCounterVar(rc));
                posIds.add(String.valueOf(rc.id()));
            }

            for(ReadContextCounter rc : phasedReadCounters.NegativeReadCounters)
            {
                negVars.add(getReadCounterVar(rc));
                negIds.add(String.valueOf(rc.id()));
            }

            SG_LOGGER.debug(String.format("LPS_DATA: %s,%s,%d,%.1f,%s,%s",
                    posIds, negIds, phasedReadCounters.ReadCount, phasedReadCounters.AllocatedReadCount, posVars, negVars));
        }
    }

    private static String getReadCounterVar(final ReadContextCounter rc)
    {
        return String.format("%d_%s>%s", rc.position(), rc.ref(), rc.alt());
    }
}
