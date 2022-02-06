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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;

public class VariantPhaser
{
    private final PhaseSetCounter mPhaseSetCounter;

    private boolean mEnabled;
    private boolean mLogData;
    private ChrBaseRegion mRegion; // only for logging

    private final List<PhasedVariantGroup> mPhasedGroups; // order by position of the lowest variant
    private int mCurrentIndex;

    public static final int MIN_READ_COUNT = 1;

    public VariantPhaser(final PhaseSetCounter phaseSetCounter)
    {
        mPhaseSetCounter = phaseSetCounter;
        mEnabled = false;
        mLogData = false;

        mPhasedGroups = Lists.newArrayList();
        mCurrentIndex = 0;
    }

    public boolean enabled() { return mEnabled; }

    public List<PhasedVariantGroup> getPhasedGroups() { return mPhasedGroups; }

    public void initialise(final ChrBaseRegion region, boolean enabled, boolean logData)
    {
        mRegion = region;
        mEnabled = enabled;
        mLogData = logData;
        mPhasedGroups.clear();
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

        /*
        if(!negCounters.isEmpty())
        {
            minVarPos = min(minVarPos, posCounters.get(0).position());
            maxVarPos = max(maxVarPos, negCounters.get(negCounters.size() - 1).position());
        }
        */

        if(mPhasedGroups.isEmpty())
        {
            mPhasedGroups.add(new PhasedVariantGroup(mPhasedGroups.size(), minVarPos, maxVarPos, posCounters, negCounters));
            return;
        }

        int index = mCurrentIndex;

        PhasedVariantGroup currentPhasedGroup = mPhasedGroups.get(index);

        //boolean searchBackwards = minVarPos < currentPhasedGroup.minVariantPos();

        boolean searchBackwards;

        if(minVarPos > currentPhasedGroup.minVariantPos())
        {
            searchBackwards = false;
        }
        else if(minVarPos < currentPhasedGroup.minVariantPos())
        {
            searchBackwards = true;
        }
        else
        {
            // may not be the first the the matching groups at this position, in which case move to the last matching
            searchBackwards = true;
            while(minVarPos == mPhasedGroups.get(index).minVariantPos())
            {
                if(index == mPhasedGroups.size() - 1)
                    break;

                index++;
            }
        }


        boolean matchFound = false;

        if(searchBackwards)
        {
            while(index >= 0)
            {
                PhasedVariantGroup phasedGroup = mPhasedGroups.get(index);

                if(phasedGroup.minVariantPos() < minVarPos)
                {
                    ++index;
                    break;
                }

                if(phasedGroup.minVariantPos() == minVarPos)
                {
                    // test for an exact match
                    if(phasedGroup.matches(minVarPos, maxVarPos, posCounters))
                    {
                        phasedGroup.ReadCount++;
                        phasedGroup.mergeNegatives(negCounters);
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
            while(index < mPhasedGroups.size())
            {
                PhasedVariantGroup phasedGroup = mPhasedGroups.get(index);

                if(minVarPos < phasedGroup.minVariantPos())
                    break;

                if(phasedGroup.minVariantPos() == minVarPos)
                {
                    // test for an exact match
                    if(phasedGroup.matches(minVarPos, maxVarPos, posCounters))
                    {
                        phasedGroup.ReadCount++;
                        phasedGroup.mergeNegatives(negCounters);
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

        if(index > 0 && index < mPhasedGroups.size() - 1)
        {
            if(minVarPos < mPhasedGroups.get(index - 1).minVariantPos() || minVarPos > mPhasedGroups.get(index + 1).minVariantPos())
            {
                SG_LOGGER.error("region({}) invalid pos() insertion at index({}) groups({})",
                        mRegion, minVarPos, index, mPhasedGroups.size());
            }
        }

        mPhasedGroups.add(index, new PhasedVariantGroup(mPhasedGroups.size(), minVarPos, maxVarPos, posCounters, negCounters));
    }

    public void mergeGroups()
    {
        /* merge groups by the following rules:
            - if a group has only matching positive read-counters (ie not a subset of any other group), then merge them all
            - if a group is a subset of 1 or more other groups, then allocate its reads pro-rata and remove it
            - finally merge any pair which share a subset (without either being a subset of each other)

        then:
            - filter for LPS support > X
            - remove any uninformative LPS (ie LPS has 1+ variant and is the only LPS that includes that variant).
        */

        /*
        List<PhasedReadCounters> phasedReadCounters = mPhasedReadCountersMap.values().stream().collect(Collectors.toList());
        Collections.sort(phasedReadCounters, new PhasedRcReadCountComparator());
        Collections.sort(phasedReadCounters, new PhasedRcPosCountComparator());
        Collections.sort(phasedReadCounters, new PhasedRcNegCountComparator());
        */

        Set<PhasedVariantGroup> removedGroups = Sets.newHashSet();

        for(PhasedVariantGroup phasedGroup : mPhasedGroups)
        {
            if(removedGroups.contains(phasedGroup))
                continue;

            // check for super-set groups
            List<PhasedVariantGroup> superGroups = mPhasedGroups.stream()
                    .filter(x -> x != phasedGroup)
                    .filter(x -> !removedGroups.contains(x))
                    .filter(x -> phasedGroup.positionsOverlap(x))
                    .filter(x -> phasedGroup.isSubsetOf(x))
                    .collect(Collectors.toList());

            if(superGroups.isEmpty())
                continue;

            // where the super groups differ, allocate by pro-rata and continue
            if(superGroups.stream().anyMatch(x -> x.PositiveReadCounters.size() > phasedGroup.PositiveReadCounters.size()))
            {
                double totalReads = superGroups.stream().mapToInt(x -> x.ReadCount).sum();

                for(PhasedVariantGroup superGroup : superGroups)
                {
                    double allocFraction = superGroup.ReadCount / totalReads;
                    superGroup.AllocatedReadCount += allocFraction * phasedGroup.ReadCount;
                    superGroup.AllocatedReadCount += allocFraction * phasedGroup.AllocatedReadCount;
                    superGroup.mergeNegatives(phasedGroup.NegativeReadCounters);
                }

                // remove subset group
                removedGroups.add(phasedGroup);
            }
            else
            {
                // otherwise collapse the groups into this one
                // this scenario is now handled during the read processing phase
                for(PhasedVariantGroup superGroup : superGroups)
                {
                    phasedGroup.ReadCount += superGroup.ReadCount;
                    phasedGroup.AllocatedReadCount += superGroup.AllocatedReadCount;
                }

                removedGroups.addAll(superGroups);
            }
        }

        removedGroups.forEach(x -> mPhasedGroups.remove(x));

        int i = 0;
        while(i < mPhasedGroups.size())
        {
            PhasedVariantGroup phasedGroup = mPhasedGroups.get(i);

            // remove any group or those with a single variant not present as a +ve in another group
            if(phasedGroup.PositiveReadCounters.size() == 1)
            {
                if(mPhasedGroups.stream()
                        .filter(x -> x != phasedGroup)
                        .filter(x -> phasedGroup.positionsOverlap(x))
                        .noneMatch(x -> x.PositiveReadCounters.contains(phasedGroup)))
                {
                    mPhasedGroups.remove(i);
                    continue;
                }
            }

            int j = i + 1;
            while(j < mPhasedGroups.size())
            {
                PhasedVariantGroup otherPhasedGroup = mPhasedGroups.get(j);

                if(!phasedGroup.positionsOverlap(otherPhasedGroup))
                    break;

                if(phasedGroup.haveCommonSubset(otherPhasedGroup))
                {
                    phasedGroup.merge(otherPhasedGroup);
                    mPhasedGroups.remove(j);
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

        int startCount = mPhasedGroups.size();

        mergeGroups();

        if(mLogData)
            logPhasedReadCounters();

        int assignedLps = 0;
        Set<ReadContextCounter> uniqueRCs = mLogData && SG_LOGGER.isTraceEnabled() ? Sets.newHashSet() : null;

        for(PhasedVariantGroup phasedGroup : mPhasedGroups)
        {
            if(phasedGroup.ReadCount + phasedGroup.AllocatedReadCount < MIN_READ_COUNT) // currently pointless since all PGs have RC >= 1
                continue;

            int nextLps = mPhaseSetCounter.getNext();
            phasedGroup.PositiveReadCounters.forEach(x -> x.addLocalPhaseSet(nextLps, phasedGroup.ReadCount, phasedGroup.AllocatedReadCount));
            ++assignedLps;

            if(uniqueRCs != null)
            {
                phasedGroup.PositiveReadCounters.forEach(x -> uniqueRCs.add(x));
                phasedGroup.NegativeReadCounters.forEach(x -> uniqueRCs.add(x));
            }
        }

        SG_LOGGER.trace("region({}) phasing groups start({}) postMerge({}) assigned({}) uniqueRCs({})",
                mRegion, startCount, mPhasedGroups.size(), assignedLps, uniqueRCs != null ? uniqueRCs.size() : 0);

        mPhasedGroups.clear();
    }

    public static class ReadCounterIdComparator implements Comparator<ReadContextCounter>
    {
        public int compare(final ReadContextCounter first, final ReadContextCounter second)
        {
            return first.id() - second.id();
        }
    }

    public static class PhasedRcReadCountComparator implements Comparator<PhasedVariantGroup>
    {
        public int compare(final PhasedVariantGroup first, final PhasedVariantGroup second)
        {
            return second.ReadCount - first.ReadCount;
        }
    }

    public static class PhasedRcPosCountComparator implements Comparator<PhasedVariantGroup>
    {
        public int compare(final PhasedVariantGroup first, final PhasedVariantGroup second)
        {
            return second.PositiveReadCounters.size() - first.PositiveReadCounters.size();
        }
    }

    public static class PhasedRcNegCountComparator implements Comparator<PhasedVariantGroup>
    {
        public int compare(final PhasedVariantGroup first, final PhasedVariantGroup second)
        {
            return second.NegativeReadCounters.size() - first.NegativeReadCounters.size();
        }
    }

    private void logPhasedReadCounters()
    {
        for(PhasedVariantGroup phasedGroup : mPhasedGroups)
        {
            StringJoiner posVars = new StringJoiner(";");
            StringJoiner posIds = new StringJoiner(";");
            StringJoiner negVars = new StringJoiner(";");
            StringJoiner negIds = new StringJoiner(";");

            for(ReadContextCounter rc : phasedGroup.PositiveReadCounters)
            {
                posVars.add(getReadCounterVar(rc));
                posIds.add(String.valueOf(rc.id()));
            }

            for(ReadContextCounter rc : phasedGroup.NegativeReadCounters)
            {
                negVars.add(getReadCounterVar(rc));
                negIds.add(String.valueOf(rc.id()));
            }

            SG_LOGGER.debug(String.format("LPS_DATA,%s,%s,%d,%.1f,%s,%s",
                    posIds, negIds, phasedGroup.ReadCount, phasedGroup.AllocatedReadCount, posVars, negVars));
        }
    }

    private static String getReadCounterVar(final ReadContextCounter rc)
    {
        return String.format("%d_%s>%s", rc.position(), rc.ref(), rc.alt());
    }
}
