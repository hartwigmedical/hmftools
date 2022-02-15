package com.hartwig.hmftools.sage.phase;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.phase.PhasedVariantGroup.maxPosition;
import static com.hartwig.hmftools.sage.phase.PhasedVariantGroup.minPosition;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public class VariantPhaser
{
    private final PhaseSetCounter mPhaseSetCounter;

    private boolean mLogData;
    private ChrBaseRegion mRegion; // only for logging

    private final List<PhasedGroupCollection> mPhasedGroupCollections; // order by position of the lowest variant
    private int mNextGroupId;

    private final List<PerformanceCounter> mPerfCounters;

    private static final int INITIAL_MIN_READ_COUNT = 1;
    private static final int FINAL_MIN_READ_COUNT = 2;
    private static final double SUBSET_READ_COUNT_LIMIT = 0.25;

    public static final int PC_PHASE_READS = 0;
    public static final int PC_FORM_LPS = 1;

    public VariantPhaser(final PhaseSetCounter phaseSetCounter)
    {
        mPhaseSetCounter = phaseSetCounter;
        mLogData = false;

        mPhasedGroupCollections = Lists.newArrayList();
        mNextGroupId = 0;

        mPerfCounters = Lists.newArrayList();
        mPerfCounters.add(new PerformanceCounter("PhaseReads"));
        mPerfCounters.add(new PerformanceCounter("FormLPS"));
    }

    public List<PhasedGroupCollection> getPhasedCollections() { return mPhasedGroupCollections; }
    public List<PerformanceCounter> getPerfCounters() { return mPerfCounters; }
    public ChrBaseRegion region() { return mRegion; }

    public void initialise(final ChrBaseRegion region, boolean logData)
    {
        mRegion = region;
        mLogData = logData;
        mPhasedGroupCollections.clear();
        mNextGroupId = 0;

        mPerfCounters.get(PC_PHASE_READS).start();
        mPerfCounters.get(PC_PHASE_READS).pause();
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
        int posVarMin = minPosition(posCounters, true);
        int posVarMax = maxPosition(posCounters, true);

        int variantMin = !negCounters.isEmpty() ? min(posVarMin, minPosition(negCounters, true)) : posVarMin;
        int variantMax = !negCounters.isEmpty() ? max(posVarMax, maxPosition(negCounters, true)) : posVarMax;

        // assign collections of variants to non-overlapping collections for subsequent comparison efficiency
        int index = 0;
        while(index < mPhasedGroupCollections.size())
        {
            PhasedGroupCollection collection = mPhasedGroupCollections.get(index);
            if(collection.positionsOverlap(variantMin, variantMax))
            {
                if(collection.addPhaseVariants(posVarMin, posVarMax, mNextGroupId, posCounters, negCounters))
                {
                    ++mNextGroupId;

                    // check for a merge with the next collection
                    if(index < mPhasedGroupCollections.size() - 1)
                    {
                        PhasedGroupCollection nextCollection = mPhasedGroupCollections.get(index + 1);
                        if(collection.positionsOverlap(nextCollection.minPosition(), nextCollection.maxPosition()))
                        {
                            collection.merge(nextCollection);
                            mPhasedGroupCollections.remove(index + 1);
                        }
                    }
                }

                return;
            }
            else if(variantMax < collection.minPosition())
            {
                break;
            }

            ++index;
        }

        PhasedGroupCollection collection = new PhasedGroupCollection(mPhasedGroupCollections.size());
        collection.addPhaseVariants(posVarMin, posVarMax, mNextGroupId, posCounters, negCounters);
        ++mNextGroupId;
        mPhasedGroupCollections.add(collection);
    }

    public void signalPhaseReadsEnd() { mPerfCounters.get(PC_PHASE_READS).stop(); }

    public void assignLocalPhaseSets(final Set<ReadContextCounter> passingCounters, final Set<ReadContextCounter> validCounters)
    {
        // assign local phase set IDs to all phased variants
        int startCount = mNextGroupId;

        boolean hasGroups = applyInitialFilters(passingCounters, validCounters);

        if(!hasGroups)
            return;

        mPerfCounters.get(PC_FORM_LPS).start();

        int startFilteredCount = mPhasedGroupCollections.stream().mapToInt(x -> x.groups().size()).sum();

        if(mLogData)
        {
            List<PhasedVariantGroup> phasedGroups = Lists.newArrayList();
            mPhasedGroupCollections.forEach(x -> phasedGroups.addAll(x.groups()));
            logPhasedReadCounters(phasedGroups, "INITIAL");
        }

        mergeGroups();

        List<PhasedVariantGroup> finalPhasedGroups = Lists.newArrayList();
        mPhasedGroupCollections.forEach(x -> finalPhasedGroups.addAll(x.groups()));

        if(mLogData)
        {
            logPhasedReadCounters(finalPhasedGroups, "FINAL");
        }

        int assignedLps = 0;
        Set<ReadContextCounter> uniqueRCs = mLogData && SG_LOGGER.isTraceEnabled() ? Sets.newHashSet() : null;

        for(PhasedVariantGroup group : finalPhasedGroups)
        {
            if(group.ReadCount + group.AllocatedReadCount < FINAL_MIN_READ_COUNT)
                continue;

            int nextLps = mPhaseSetCounter.getNext();
            group.PositiveReadCounters.forEach(x -> x.addLocalPhaseSet(nextLps, group.ReadCount, group.AllocatedReadCount));
            ++assignedLps;

            if(uniqueRCs != null)
            {
                group.PositiveReadCounters.forEach(x -> uniqueRCs.add(x));
                group.NegativeReadCounters.forEach(x -> uniqueRCs.add(x));
            }
        }

        SG_LOGGER.trace("region({}) phasing groups(coll={} start={} filtered={}) postMerge({}) assigned({}) rc(pass={} valid={} uniqueRCs={})",
                mRegion, mPhasedGroupCollections.size(), startCount, startFilteredCount, finalPhasedGroups.size(), assignedLps,
                passingCounters.size(), validCounters.size(), uniqueRCs != null ? uniqueRCs.size() : 0);

        mPerfCounters.get(PC_FORM_LPS).stop();

        if(SG_LOGGER.isDebugEnabled())
        {
            PerformanceCounter pc = mPerfCounters.get(PC_FORM_LPS);
            double lastTime = pc.getTimes().get(pc.getTimes().size() - 1);
            if(lastTime > 0.5)
            {
                SG_LOGGER.debug("region({}) phasing groups start({} filtered={}) postMerge({}) assigned({}) rc(pass={} valid={}) time({})",
                        mRegion, startCount, startFilteredCount, finalPhasedGroups.size(), assignedLps,
                        passingCounters.size(), validCounters.size(), lastTime);
            }
        }
    }

    private boolean applyInitialFilters(
            final Set<ReadContextCounter> passingReadCounters, final Set<ReadContextCounter> validCounters)
    {
        // returns true if a group has variants left after filtering
        if(passingReadCounters.isEmpty())
            return false;

        boolean hasGroups = false;

        for(PhasedGroupCollection collection : mPhasedGroupCollections)
        {
            int index = 0;
            while(index < collection.groups().size())
            {
                PhasedVariantGroup group = collection.groups().get(index);

                if(group.ReadCount < INITIAL_MIN_READ_COUNT)
                {
                    collection.groups().remove(index);
                    continue;
                }

                if(group.PositiveReadCounters.stream().noneMatch(y -> passingReadCounters.contains(y)))
                {
                    collection.groups().remove(index);
                    continue;
                }

                if(group.cullReadCounters(validCounters) && !group.isValid())
                {
                    collection.groups().remove(index);
                    continue;
                }

                ++index;
            }

            hasGroups |= !collection.groups().isEmpty();
        }

        return hasGroups;
    }

    public void mergeGroups()
    {
        // reduce groups to minimum set of information local phase sets by merging:
        // groups with matching positives, non-conflicting negatives and not subsets of other groups
        // extending groups with common overlaps if only one group extends either up or down
        // removing uninformative groups

        for(PhasedGroupCollection collection : mPhasedGroupCollections)
        {
            // then apply merging rules within these overlapping groups
            mergeMatching(collection.groups());

            mergeByExtension(collection.groups());

            mergeUninformative(collection.groups());

            mergeByExtension(collection.groups());
        }
    }

    public static void mergeMatching(final List<PhasedVariantGroup> filteredGroups)
    {
        // merge groups with matching positives, non-conflicting negatives and not subsets of other groups
        if(filteredGroups.size() < 2)
            return;

        int i = 0;
        while(i < filteredGroups.size())
        {
            PhasedVariantGroup group = filteredGroups.get(i);

            // find groups which have matching +ves and aren't subsets of any other group
            List<PhasedVariantGroup> superGroups = Lists.newArrayList();

            for(int j = i + 1; j < filteredGroups.size(); ++j)
            {
                PhasedVariantGroup otherGroup = filteredGroups.get(j);

                if(group.positionsOverlap(otherGroup) && group.isSubsetOf(otherGroup))
                    superGroups.add(otherGroup);
            }

            if(superGroups.isEmpty())
            {
                if(group.PositiveReadCounters.size() == 1)
                {
                    // remove any group with a single variant not present as a +ve in another group
                    filteredGroups.remove(group);
                    continue;
                }
            }
            else
            {
                // if all groups have matching +ves then collapse into the current
                if(superGroups.stream().noneMatch(x -> x.PositiveReadCounters.size() > group.PositiveReadCounters.size()))
                {
                    for(PhasedVariantGroup superGroup : superGroups)
                    {
                        group.merge(superGroup);
                    }

                    filteredGroups.removeAll(superGroups);
                }
            }

            ++i;
        }
    }

    public static void mergeByExtension(final List<PhasedVariantGroup> filteredGroups)
    {
        if(filteredGroups.size() < 2)
            return;

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
                PhasedVariantGroup group = filteredGroups.get(i);

                if(!initialLoop && !lastModifiedGroups.contains(group))
                    continue;

                List<PhasedVariantGroup> candidateGroups = Lists.newArrayList(group);

                commonPosCounters.clear();
                commonNegCounters.clear();

                for(int j = 0; j < filteredGroups.size(); ++j)
                {
                    if(j == i)
                        continue;

                    PhasedVariantGroup otherPhasedGroup = filteredGroups.get(j);

                    if(!group.positionsOverlap(otherPhasedGroup))
                        continue;

                    if(commonPosCounters.isEmpty())
                    {
                        if(!group.populateCommon(otherPhasedGroup, commonPosCounters, commonNegCounters))
                            continue;

                        candidateGroups.add(otherPhasedGroup);
                    }
                    //else if(otherPhasedGroup.hasCommonSubset(group, commonPosCounters, commonNegCounters))
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

    public static void mergeUninformative(final List<PhasedVariantGroup> filteredGroups)
    {
        // finally merge any groups with the same +ves or are non-conflicting subsets of others now that supersets have been considered
        int index = 0;
        while(index < filteredGroups.size())
        {
            PhasedVariantGroup group = filteredGroups.get(index);

            if(group.PositiveReadCounters.size() == 1)
            {
                ReadContextCounter readCounter = group.PositiveReadCounters.get(0);

                if(filteredGroups.stream().filter(x -> x != group).noneMatch(x -> group.PositiveReadCounters.contains(readCounter)))
                {
                    // remove any group with a single variant not present as a +ve in another group
                    filteredGroups.remove(group);
                    continue;
                }
            }

            // find groups which have matching +ves and aren't subsets of any other group
            List<PhasedVariantGroup> matchingGroups = filteredGroups.stream()
                    .filter(x -> x != group)
                    .filter(x -> group.positivesMatch(x) || x.isSubsetOf(group))
                    .collect(Collectors.toList());

            if(!matchingGroups.isEmpty())
            {
                // collapse the groups into this one
                for(PhasedVariantGroup otherGroup : matchingGroups)
                {
                    group.merge(otherGroup);
                }

                filteredGroups.removeAll(matchingGroups);
            }

            ++index;
        }
    }

    public static void removeUninformativeLps(final List<SageVariant> variants, final Set<Integer> passingPhaseSets)
    {
        // remove any uninformative local phasings sets where they all have the same passing variants
        Map<Integer,List<SageVariant>> lpsVariantsMap = Maps.newHashMap();
        Map<Integer,List<SageVariant>> lpsPassingVariantsMap = Maps.newHashMap();
        Map<Integer,Integer> lpsMaxReadCountMap = Maps.newHashMap();

        List<Integer> uninformativeLpsIds = Lists.newArrayList();
        List<Integer> processedLpsIds = Lists.newArrayList();
        Set<Integer> singlePassingVarGroups = Sets.newHashSet();

        // first put all variants into LPS datasets
        for(Integer lpsId : passingPhaseSets)
        {
            List<SageVariant> lpsVariants = variants.stream().filter(x -> x.hasMatchingLps(lpsId)).collect(Collectors.toList());
            List<SageVariant> passingVariants = lpsVariants.stream().filter(x -> x.isPassing()).collect(Collectors.toList());

            if(passingVariants.isEmpty())
            {
                uninformativeLpsIds.add(lpsId);
                processedLpsIds.add(lpsId);
                continue;
            }

            if(lpsVariants.size() == 1 && passingVariants.get(0).localPhaseSets().size() == 1)
                singlePassingVarGroups.add(lpsId);

            lpsVariantsMap.put(lpsId, lpsVariants);
            lpsPassingVariantsMap.put(lpsId, passingVariants);
            lpsMaxReadCountMap.put(lpsId, lpsVariants.get(0).getLpsReadCount(lpsId));
        }

        for(Map.Entry<Integer,List<SageVariant>> entry : lpsVariantsMap.entrySet())
        {
            Integer lpsId = entry.getKey();

            if(processedLpsIds.contains(lpsId))
                continue;

            processedLpsIds.add(lpsId);

            List<SageVariant> passingVariants = lpsPassingVariantsMap.get(lpsId);
            int maxReadCount = passingVariants.get(0).getLpsReadCount(lpsId);

            // look for matching passing variants
            int maxLpsId = lpsId;

            List<Integer> matchedVariantsLpsIds = Lists.newArrayList(lpsId);

            for(Map.Entry<Integer,List<SageVariant>> entry2 : lpsVariantsMap.entrySet())
            {
                Integer otherLpsId = entry2.getKey();

                if(otherLpsId == lpsId || processedLpsIds.contains(otherLpsId))
                    continue;

                // List<SageVariant> otherLpsVariants = entry.getValue();
                List<SageVariant> otherPassingVariants = lpsPassingVariantsMap.get(otherLpsId);

                if(otherPassingVariants.size() == passingVariants.size() && otherPassingVariants.stream().allMatch(x -> passingVariants.contains(x)))
                {
                    processedLpsIds.add(otherLpsId);

                    int otherReadCount = otherPassingVariants.get(0).getLpsReadCount(otherLpsId);
                    matchedVariantsLpsIds.add(otherLpsId);

                    if(otherReadCount > maxReadCount)
                    {
                        maxReadCount = otherReadCount;
                        maxLpsId = otherLpsId;
                    }
                }
            }

            if(matchedVariantsLpsIds.size() > 1)
            {
                int maxId = maxLpsId;
                matchedVariantsLpsIds.stream().filter(x -> x != maxId).forEach(x -> uninformativeLpsIds.add(x));

                if(passingVariants.size() == 1)
                    singlePassingVarGroups.add(maxLpsId);
            }
        }

        // and filter any subsets of PASS variants with <25% read support of a superset of PASS variants
        for(Map.Entry<Integer,List<SageVariant>> entry : lpsVariantsMap.entrySet())
        {
            Integer lpsId = entry.getKey();

            if(uninformativeLpsIds.contains(lpsId))
                continue;

            List<SageVariant> passingVariants = lpsPassingVariantsMap.get(lpsId);

            int readCount = passingVariants.get(0).getLpsReadCount(lpsId);

            // look for superset of these passing variants
            for(Map.Entry<Integer, List<SageVariant>> entry2 : lpsVariantsMap.entrySet())
            {
                Integer otherLpsId = entry2.getKey();

                if(otherLpsId == lpsId || uninformativeLpsIds.contains(otherLpsId))
                    continue;

                List<SageVariant> otherPassingVariants = lpsPassingVariantsMap.get(otherLpsId);

                if(otherPassingVariants.size() > passingVariants.size() && passingVariants.stream().allMatch(x -> otherPassingVariants.contains(x)))
                {
                    // is a subset
                    int otherReadCount = otherPassingVariants.get(0).getLpsReadCount(otherLpsId);
                    if(readCount < otherReadCount * SUBSET_READ_COUNT_LIMIT)
                    {
                        uninformativeLpsIds.add(lpsId);

                        if(otherPassingVariants.size() == 1)
                            singlePassingVarGroups.add(otherLpsId);

                        break;
                    }
                }
            }
        }

        // now remove all uninformative IDs from each variant
        for(Integer lpsId : uninformativeLpsIds)
        {
            List<SageVariant> lpsVariants = lpsVariantsMap.get(lpsId);

            lpsVariants.forEach(x -> x.removeLps(lpsId));
        }

        // remove any LPS if a variant now has no other LPS and is by itself
        for(Integer lpsId : singlePassingVarGroups)
        {
            List<SageVariant> lpsVariants = lpsVariantsMap.get(lpsId);

            if(lpsVariants.stream().allMatch(x -> x.hasLocalPhaseSets() && x.localPhaseSets().size() == 1))
            {
                lpsVariants.forEach(x -> x.removeLps(lpsId));
            }
        }
    }

    private void logPhasedReadCounters(final List<PhasedVariantGroup> phasedGroups, final String stage)
    {
        for(PhasedVariantGroup group : phasedGroups)
        {
            StringJoiner posVars = new StringJoiner(";");
            StringJoiner posIds = new StringJoiner(";");
            StringJoiner negVars = new StringJoiner(";");
            StringJoiner negIds = new StringJoiner(";");

            for(ReadContextCounter rc : group.PositiveReadCounters)
            {
                posVars.add(getReadCounterVar(rc));
                posIds.add(String.valueOf(rc.id()));
            }

            for(ReadContextCounter rc : group.NegativeReadCounters)
            {
                negVars.add(getReadCounterVar(rc));
                negIds.add(String.valueOf(rc.id()));
            }

            // Time,Stage,Id,MergedIds,PosIds,NegIds,UniqueReadCount,AllocReadCount,PosVars,NegVars
            SG_LOGGER.debug(String.format("LPS_DATA,%s,%d,%s,%s,%s,%d,%.1f,%s,%s",
                    stage, group.Id, group.mergedGroupIds(), posIds, negIds,
                    group.ReadCount, group.AllocatedReadCount, posVars, negVars));
        }
    }

    private static String getReadCounterVar(final ReadContextCounter rc)
    {
        return String.format("%d_%s>%s", rc.position(), rc.ref(), rc.alt());
    }

}
