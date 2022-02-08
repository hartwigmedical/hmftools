package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
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

    public static final int INITIAL_MIN_READ_COUNT = 1;
    public static final int FINAL_MIN_READ_COUNT = 2;

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

        if(minVarPos > currentPhasedGroup.posVariantMin())
        {
            searchBackwards = false;
        }
        else if(minVarPos < currentPhasedGroup.posVariantMin())
        {
            searchBackwards = true;
        }
        else
        {
            // may not be the first the the matching groups at this position, in which case move to the last matching
            searchBackwards = true;
            while(minVarPos == mPhasedGroups.get(index).posVariantMin())
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

                if(phasedGroup.posVariantMin() < minVarPos)
                {
                    ++index;
                    break;
                }

                if(phasedGroup.posVariantMin() == minVarPos)
                {
                    // test for an exact match
                    if(phasedGroup.exactMatch(minVarPos, maxVarPos, posCounters, negCounters))
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

                if(minVarPos < phasedGroup.posVariantMin())
                    break;

                if(phasedGroup.posVariantMin() == minVarPos)
                {
                    // test for an exact match
                    if(phasedGroup.exactMatch(minVarPos, maxVarPos, posCounters, negCounters))
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
            if(minVarPos < mPhasedGroups.get(index - 1).posVariantMin() || minVarPos > mPhasedGroups.get(index + 1).posVariantMin())
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
            -

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

        List<PhasedVariantGroup> filteredGroups = mPhasedGroups.stream()
                .filter(x -> x.ReadCount >= INITIAL_MIN_READ_COUNT).collect(Collectors.toList());

        mergeMatching(filteredGroups);

        mergeByExtension(filteredGroups);

        mergeUninformative(filteredGroups);

        mergeByExtension(filteredGroups);

        mPhasedGroups.clear();
        mPhasedGroups.addAll(filteredGroups);
    }

    private void mergeMatching(final List<PhasedVariantGroup> filteredGroups)
    {
        Set<PhasedVariantGroup> removedGroups = Sets.newHashSet();

        for(PhasedVariantGroup phasedGroup : filteredGroups)
        {
            if(removedGroups.contains(phasedGroup))
                continue;

            // find groups which have matching +ves and aren't subsets of any other group
            List<PhasedVariantGroup> superGroups = filteredGroups.stream()
                    .filter(x -> x != phasedGroup)
                    .filter(x -> !removedGroups.contains(x))
                    .filter(x -> phasedGroup.positionsOverlap(x))
                    .filter(x -> phasedGroup.isSubsetOf(x))
                    .collect(Collectors.toList());

            if(superGroups.isEmpty())
            {
                // remove any group with a single variant not present as a +ve in another group
                if(phasedGroup.PositiveReadCounters.size() == 1)
                    removedGroups.add(phasedGroup);

                continue;
            }

            if(superGroups.stream().noneMatch(x -> x.PositiveReadCounters.size() > phasedGroup.PositiveReadCounters.size()))
            {
                // collapse the groups into this one
                for(PhasedVariantGroup superGroup : superGroups)
                {
                    phasedGroup.merge(superGroup);
                }

                removedGroups.addAll(superGroups);
            }
        }

        removedGroups.forEach(x -> filteredGroups.remove(x));
    }

    private void mergeByExtension(final List<PhasedVariantGroup> filteredGroups)
    {
        // merge any group with a common subset of +ves and -ves if it can only be extended in one direction
        List<ReadContextCounter> commonPosCounters = Lists.newArrayList();
        List<ReadContextCounter> commonNegCounters = Lists.newArrayList();
        Set<PhasedVariantGroup> modifiedGroups = Sets.newHashSet();
        Set<PhasedVariantGroup> lastModifiedGroups = Sets.newHashSet();
        boolean initialLoop = true;

        while(initialLoop || !modifiedGroups.isEmpty())
        {
            lastModifiedGroups.clear();
            lastModifiedGroups.addAll(modifiedGroups);
            modifiedGroups.clear();

            for(int i = 0; i < filteredGroups.size(); ++i)
            {
                PhasedVariantGroup phasedGroup = filteredGroups.get(i);

                if(!initialLoop && !lastModifiedGroups.contains(phasedGroup))
                    continue;

                List<PhasedVariantGroup> candidateGroups = Lists.newArrayList(phasedGroup);

                commonPosCounters.clear();
                commonNegCounters.clear();

                for(int j = i + 1; j < filteredGroups.size(); ++j)
                {
                    PhasedVariantGroup otherPhasedGroup = filteredGroups.get(j);

                    if(!phasedGroup.positionsOverlap(otherPhasedGroup))
                        continue;

                    if(commonPosCounters.isEmpty())
                    {
                        if(!phasedGroup.populateCommon(otherPhasedGroup, commonPosCounters, commonNegCounters))
                            continue;

                        candidateGroups.add(otherPhasedGroup);
                    }
                    //else if(otherPhasedGroup.hasCommonSubset(phasedGroup, commonPosCounters, commonNegCounters))
                    else if(candidateGroups.stream().allMatch(x -> otherPhasedGroup.hasCommonSubset(x, commonPosCounters, commonNegCounters)))
                    {
                        candidateGroups.add(otherPhasedGroup);
                    }
                }

                if(candidateGroups.size() == 1)
                    continue;

                // check for a single group that extends the common subset in one direction
                commonPosCounters.addAll(commonNegCounters);
                int minSubsetPos = commonPosCounters.stream().mapToInt(x -> x.position()).min().orElse(0);
                int maxSubsetPos = commonPosCounters.stream().mapToInt(x -> x.position()).max().orElse(0);

                List<PhasedVariantGroup> lowerGroups = candidateGroups.stream()
                        .filter(x -> x.variantMin() < minSubsetPos && x.variantMax() == maxSubsetPos).collect(Collectors.toList());

                List<PhasedVariantGroup> upperGroups = candidateGroups.stream()
                        .filter(x -> x.variantMin() == minSubsetPos && x.variantMax() > maxSubsetPos).collect(Collectors.toList());

                if(lowerGroups.size() == 1 || upperGroups.size() == 1)
                {
                    PhasedVariantGroup mergedGroup = lowerGroups.size() == 1 ? lowerGroups.get(0) : upperGroups.get(0);

                    List<PhasedVariantGroup> mergingGroups = lowerGroups.size() == 1 ? upperGroups : lowerGroups;

                    if(mergingGroups.isEmpty())
                        mergingGroups = candidateGroups.stream().filter(x -> x != mergedGroup).collect(Collectors.toList());

                    // merge into others
                    double totalReads = mergingGroups.stream().mapToInt(x -> x.ReadCount).sum();

                    for(PhasedVariantGroup otherGroup : mergingGroups)
                    {
                        double allocFraction = otherGroup.ReadCount / totalReads;
                        otherGroup.merge(mergedGroup, allocFraction);
                        modifiedGroups.add(otherGroup);
                    }

                    filteredGroups.remove(mergedGroup);
                    break;
                }
            }

            initialLoop = false;
        }
    }

    private void mergeUninformative(final List<PhasedVariantGroup> filteredGroups)
    {
        // finally merge any groups with the same +ves or are non-conflicting subsets of others now that supersets have been considered
        int index = 0;
        while(index < filteredGroups.size())
        {
            PhasedVariantGroup phasedGroup = filteredGroups.get(index);

            // find groups which have matching +ves and aren't subsets of any other group
            List<PhasedVariantGroup> matchingGroups = filteredGroups.stream()
                    .filter(x -> x != phasedGroup)
                    .filter(x -> phasedGroup.positivesMatch(x) || x.isSubsetOf(phasedGroup))
                    .collect(Collectors.toList());

            if(!matchingGroups.isEmpty())
            {
                // collapse the groups into this one
                for(PhasedVariantGroup otherGroup : matchingGroups)
                {
                    phasedGroup.merge(otherGroup);
                }

                filteredGroups.removeAll(matchingGroups);
            }

            ++index;
        }
    }

    public void assignLocalPhaseSets()
    {
        // assign local phase set IDs to all phased variants
        if(!mEnabled)
            return;

        int startCount = mPhasedGroups.size();

        if(mLogData)
            logPhasedReadCounters(mPhasedGroups, "INITIAL");

        mergeGroups();

        if(mLogData)
            logPhasedReadCounters(mPhasedGroups, "FINAL");

        int assignedLps = 0;
        Set<ReadContextCounter> uniqueRCs = mLogData && SG_LOGGER.isTraceEnabled() ? Sets.newHashSet() : null;

        for(PhasedVariantGroup phasedGroup : mPhasedGroups)
        {
            if(phasedGroup.ReadCount + phasedGroup.AllocatedReadCount < FINAL_MIN_READ_COUNT)
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

    private void logPhasedReadCounters(final List<PhasedVariantGroup> phasedGroups, final String stage)
    {
        for(PhasedVariantGroup phasedGroup : phasedGroups)
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

            // Time,Stage,Id,MergedIds,PosIds,NegIds,UniqueReadCount,AllocReadCount,PosVars,NegVars
            SG_LOGGER.debug(String.format("LPS_DATA,%s,%d,%s,%s,%s,%d,%.1f,%s,%s",
                    stage, phasedGroup.Id, phasedGroup.mergedGroupIds(), posIds, negIds,
                    phasedGroup.ReadCount, phasedGroup.AllocatedReadCount, posVars, negVars));
        }
    }

    private static String getReadCounterVar(final ReadContextCounter rc)
    {
        return String.format("%d_%s>%s", rc.position(), rc.ref(), rc.alt());
    }
}
