package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;

public class VariantPhaser
{
    private final PhaseSetCounter mPhaseSetCounter;

    private boolean mLogData;
    private ChrBaseRegion mRegion; // only for logging

    private final List<PhasedVariantGroup> mPhasedGroups; // order by position of the lowest variant
    private int mCurrentIndex;

    private final List<PerformanceCounter> mPerfCounters;

    public static final int INITIAL_MIN_READ_COUNT = 1;
    public static final int FINAL_MIN_READ_COUNT = 2;
    public static final int PC_PHASE_READS = 0;
    public static final int PC_FORM_LPS = 1;

    public VariantPhaser(final PhaseSetCounter phaseSetCounter)
    {
        mPhaseSetCounter = phaseSetCounter;
        mLogData = false;

        mPhasedGroups = Lists.newArrayList();
        mCurrentIndex = 0;

        mPerfCounters = Lists.newArrayList();
        mPerfCounters.add(new PerformanceCounter("PhaseReads"));
        mPerfCounters.add(new PerformanceCounter("FormLPS"));
    }

    public List<PhasedVariantGroup> getPhasedGroups() { return mPhasedGroups; }
    public List<PerformanceCounter> getPerfCounters() { return mPerfCounters; }

    public void initialise(final ChrBaseRegion region, boolean logData)
    {
        mRegion = region;
        mLogData = logData;
        mPhasedGroups.clear();
        mCurrentIndex = 0;

        mPerfCounters.get(PC_PHASE_READS).start();
        mPerfCounters.get(PC_PHASE_READS).resume();
    }

    public void registeredPhasedVariants(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        if(posCounters.isEmpty() || posCounters.size() + negCounters.size() < 2)
            return;

        mPerfCounters.get(PC_PHASE_READS).resume();
        processPhasedVariants(posCounters, negCounters);
        mPerfCounters.get(PC_PHASE_READS).pause();
    }

    private void processPhasedVariants(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        // walk backwards then forwards from the current location looking for a match
        int posVarMin = posCounters.get(0).position();
        int posVarMax = posCounters.get(posCounters.size() - 1).position();

        if(mPhasedGroups.isEmpty())
        {
            mPhasedGroups.add(new PhasedVariantGroup(mPhasedGroups.size(), posVarMin, posVarMax, posCounters, negCounters));
            return;
        }

        int index = mCurrentIndex;

        PhasedVariantGroup currentPhasedGroup = mPhasedGroups.get(index);

        //boolean searchBackwards = minVarPos < currentPhasedGroup.minVariantPos();

        boolean searchBackwards;

        if(posVarMin > currentPhasedGroup.posVariantMin())
        {
            searchBackwards = false;
        }
        else if(posVarMin < currentPhasedGroup.posVariantMin())
        {
            searchBackwards = true;
        }
        else
        {
            // may not be the first the the matching groups at this position, in which case move to the last matching
            searchBackwards = true;
            while(posVarMin == mPhasedGroups.get(index).posVariantMin())
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

                if(phasedGroup.posVariantMin() < posVarMin)
                {
                    ++index;
                    break;
                }

                if(phasedGroup.posVariantMin() == posVarMin)
                {
                    // test for an exact match
                    if(phasedGroup.exactMatch(posVarMin, posVarMax, posCounters, negCounters))
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

                if(posVarMin < phasedGroup.posVariantMin())
                    break;

                if(phasedGroup.posVariantMin() == posVarMin)
                {
                    // test for an exact match
                    if(phasedGroup.exactMatch(posVarMin, posVarMax, posCounters, negCounters))
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
            if(posVarMin < mPhasedGroups.get(index - 1).posVariantMin() || posVarMin > mPhasedGroups.get(index + 1).posVariantMin())
            {
                SG_LOGGER.error("region({}) invalid pos() insertion at index({}) groups({})",
                        mRegion, posVarMin, index, mPhasedGroups.size());
            }
        }

        mPhasedGroups.add(index, new PhasedVariantGroup(mPhasedGroups.size(), posVarMin, posVarMax, posCounters, negCounters));
    }

    public void assignLocalPhaseSets(final Set<ReadContextCounter> passingReadCounters)
    {
        // assign local phase set IDs to all phased variants
        mPerfCounters.get(PC_PHASE_READS).stop();

        int startCount = mPhasedGroups.size();

        List<PhasedVariantGroup> filteredGroups = applyInitialFilters(passingReadCounters);

        mPhasedGroups.clear();

        if(filteredGroups.isEmpty())
            return;

        mPerfCounters.get(PC_FORM_LPS).start();

        int startFilteredCount = filteredGroups.size();

        if(mLogData)
            logPhasedReadCounters(filteredGroups, "INITIAL");

        mergeGroups(filteredGroups);

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

        SG_LOGGER.trace("region({}) phasing groups start({} filtered={}) postMerge({}) assigned({}) uniqueRCs({})",
                mRegion, startCount, startFilteredCount, mPhasedGroups.size(), assignedLps, uniqueRCs != null ? uniqueRCs.size() : 0);

        mPhasedGroups.clear();

        mPerfCounters.get(PC_FORM_LPS).stop();
    }

    private List<PhasedVariantGroup> applyInitialFilters(final Set<ReadContextCounter> passingReadCounters)
    {
        List<PhasedVariantGroup> filteredGroups = mPhasedGroups.stream()
                .filter(x -> x.ReadCount >= INITIAL_MIN_READ_COUNT)
                .filter(x -> x.PositiveReadCounters.stream().anyMatch(y -> passingReadCounters.contains(y)))
                .collect(Collectors.toList());

        return filteredGroups;
    }

    public void mergeGroups(final List<PhasedVariantGroup> filteredGroups)
    {
        /* merge groups by the following rules:
            - if a group has only matching positive read-counters (ie not a subset of any other group), then merge them all
            -

        then:
            - filter for LPS support > X
            - remove any uninformative LPS (ie LPS has 1+ variant and is the only LPS that includes that variant).
        */

        // break groups into non-overlapping collections
        while(!filteredGroups.isEmpty())
        {
            List<PhasedVariantGroup> overlappingGroups = Lists.newArrayList(filteredGroups.get(0));
            int groupMin = filteredGroups.get(0).variantMin();
            int groupMax = filteredGroups.get(0).variantMax();

            filteredGroups.remove(0);

            int index = 0;
            while(index < filteredGroups.size())
            {
                PhasedVariantGroup phasedGroup = filteredGroups.get(index);

                if(positionsOverlap(groupMin, groupMax, phasedGroup.variantMin(), phasedGroup.variantMax()))
                {
                    groupMin = min(groupMin, phasedGroup.variantMin());
                    groupMax = max(groupMax, phasedGroup.variantMax());
                    overlappingGroups.add(phasedGroup);
                    filteredGroups.remove(index);
                }
                else
                {
                    ++index;
                }
            }

            // then apply merging rules within these overlapping groups
            mergeMatching(overlappingGroups);

            mergeByExtension(overlappingGroups);

            mergeUninformative(overlappingGroups);

            mergeByExtension(overlappingGroups);

            mPhasedGroups.addAll(overlappingGroups);
        }
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
