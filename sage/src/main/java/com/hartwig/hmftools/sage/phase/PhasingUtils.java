package com.hartwig.hmftools.sage.phase;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public final class PhasingUtils
{
    private static final double SUBSET_READ_COUNT_LIMIT = 0.25;

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

            lpsVariantsMap.put(lpsId, lpsVariants);

            if(passingVariants.isEmpty())
            {
                uninformativeLpsIds.add(lpsId);
                processedLpsIds.add(lpsId);
                continue;
            }

            if(lpsVariants.size() == 1 && passingVariants.get(0).localPhaseSets().size() == 1)
                singlePassingVarGroups.add(lpsId);

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
}
