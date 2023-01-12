package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.calcFragmentStatus;
import static com.hartwig.hmftools.bamtools.markdups.UmiGroup.exceedsUmiIdDiff;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroupUtils
{
    private static final int HIGH_DEPTH_THRESHOLD = 10000;

    public static void findDuplicateFragments(
            final List<Fragment> fragments, final List<Fragment> resolvedFragments,
            final List<List<Fragment>> duplicateGroups, final List<CandidateDuplicates> candidateDuplicatesList)
    {
        // take all the fragments at this initial fragment position and classify them as duplicates, non-duplicates (NONE) or unclear
        // note: all fragments will be given a classification, and resolved fragments are removed from the input fragment list

        if(fragments.size() == 1)
        {
            Fragment fragment = fragments.get(0);
            fragment.setStatus(NONE);
            resolvedFragments.add(fragment);
            fragments.clear();
            return;
        }

        int fragmentCount = fragments.size();
        boolean applyHighDepthLogic = false; // fragmentCount > HIGH_DEPTH_THRESHOLD;

        Map<Fragment,List<Fragment>> possibleDuplicates = Maps.newHashMap();
        Map<Fragment,Fragment> linkedDuplicates = Maps.newHashMap();

        int i = 0;
        while(i < fragments.size())
        {
            Fragment fragment1 = fragments.get(i);

            if(i == fragments.size() - 1)
            {
                if(!possibleDuplicates.containsKey(fragment1) && !linkedDuplicates.containsKey(fragment1))
                {
                    fragment1.setStatus(NONE);
                    resolvedFragments.add(fragment1);
                    fragments.remove(i);
                }
                break;
            }

            List<Fragment> duplicateFragments = null;

            Fragment existingLinkedFragment1 = linkedDuplicates.get(fragment1);

            boolean isCandidateDup = existingLinkedFragment1 != null;

            List<Fragment> candidateFragments = existingLinkedFragment1 != null ? possibleDuplicates.get(existingLinkedFragment1) : null;

            int j = i + 1;
            while(j < fragments.size())
            {
                Fragment fragment2 = fragments.get(j);

                FragmentStatus status = calcFragmentStatus(fragment1, fragment2);

                if(applyHighDepthLogic && status == CANDIDATE && proximateFragmentSizes(fragment1, fragment2))
                    status = DUPLICATE;

                if(status == DUPLICATE)
                {
                    fragment1.setStatus(status);
                    fragment2.setStatus(status);

                    if(duplicateFragments == null)
                        duplicateFragments = Lists.newArrayList(fragment1);

                    duplicateFragments.add(fragment2);
                    fragments.remove(j);
                    continue;
                }

                if(fragment1.status() != DUPLICATE && status == CANDIDATE)
                {
                    isCandidateDup = true;

                    // the pair is a candidate for duplicates but without their mates it's unclear whether they will be
                    Fragment existingLinkedFragment2 = linkedDuplicates.get(fragment2);

                    if(existingLinkedFragment1 != null && existingLinkedFragment1 == existingLinkedFragment2)
                    {
                        // already a part of the same candidate group
                    }
                    else if(existingLinkedFragment2 != null)
                    {
                        List<Fragment> existingGroup = possibleDuplicates.get(existingLinkedFragment2);

                        if(candidateFragments == null)
                        {
                            existingGroup.add(fragment1);
                        }
                        else
                        {
                            // take this fragment's candidates and move them to the existing group
                            for(Fragment fragment : candidateFragments)
                            {
                                if(!existingGroup.contains(fragment))
                                    existingGroup.add(fragment);

                                linkedDuplicates.put(fragment, existingLinkedFragment2);
                            }

                            possibleDuplicates.remove(existingLinkedFragment1);
                        }

                        linkedDuplicates.put(fragment1, existingLinkedFragment2);
                        existingLinkedFragment1 = existingLinkedFragment2;
                        candidateFragments = existingGroup;
                    }
                    else
                    {
                        if(candidateFragments == null)
                        {
                            candidateFragments = Lists.newArrayList(fragment1);
                            possibleDuplicates.put(fragment1, candidateFragments);
                            existingLinkedFragment1 = fragment1;
                        }

                        candidateFragments.add(fragment2);
                        linkedDuplicates.put(fragment2, existingLinkedFragment1);
                    }
                }

                ++j;
            }

            if(fragment1.status().isDuplicate())
            {
                if(isCandidateDup && possibleDuplicates.containsKey(fragment1))
                {
                    // clean-up
                    candidateFragments.forEach(x -> linkedDuplicates.remove(x));
                    possibleDuplicates.remove(fragment1);
                }

                int dupCount = duplicateFragments.size();
                duplicateFragments.forEach(x -> x.setDuplicateCount(dupCount));

                resolvedFragments.addAll(duplicateFragments);
                fragments.remove(i);

                duplicateGroups.add(duplicateFragments);
            }
            else if(isCandidateDup)
            {
                ++i;
            }
            else
            {
                fragment1.setStatus(NONE);
                resolvedFragments.add(fragment1);
                fragments.remove(i);
            }
        }

        if(!possibleDuplicates.isEmpty())
        {
            for(List<Fragment> candidateFragments : possibleDuplicates.values())
            {
                List<Fragment> unresolvedFragments = candidateFragments.stream().filter(x -> !x.status().isResolved()).collect(Collectors.toList());

                if(unresolvedFragments.size() >= 2)
                {
                    CandidateDuplicates candidateDuplicates = CandidateDuplicates.from(unresolvedFragments.get(0));
                    candidateDuplicatesList.add(candidateDuplicates);

                    for(int index = 1; index < unresolvedFragments.size(); ++index)
                    {
                        candidateDuplicates.addFragment(unresolvedFragments.get(index));
                    }
                }
                else if(unresolvedFragments.size() == 1)
                {
                    Fragment fragment = unresolvedFragments.get(0);
                    fragment.setStatus(NONE);
                    resolvedFragments.add(fragment);
                }
            }
        }

        int candidateCount = 0;
        for(CandidateDuplicates candidateDuplicates : candidateDuplicatesList)
        {
            candidateDuplicates.fragments().forEach(x -> x.setStatus(CANDIDATE));
            candidateDuplicates.fragments().forEach(x -> x.setCandidateDupKey(candidateDuplicates.key()));
            candidateCount += candidateDuplicates.fragmentCount();
        }

        if(candidateCount + resolvedFragments.size() != fragmentCount)
        {
            BM_LOGGER.error("failed to classify all fragments: original({}) resolved({}) candidates({})",
                    fragmentCount, resolvedFragments.size(), candidateCount);

            for(Fragment fragment : resolvedFragments)
            {
                BM_LOGGER.error("resolved fragment: {}", fragment);
            }

            for(CandidateDuplicates candidateDuplicates : candidateDuplicatesList)
            {
                for(Fragment fragment : candidateDuplicates.fragments())
                {
                    BM_LOGGER.error("candidate dup fragment: {}", fragment);
                }
            }
        }
    }

    private static boolean proximateFragmentSizes(final Fragment first, final Fragment second)
    {
        return abs(abs(first.reads().get(0).getInferredInsertSize()) - abs(second.reads().get(0).getInferredInsertSize())) <= 10;
    }

    public static void processDuplicateGroups(
            final List<List<Fragment>> duplicateGroups, final UmiConfig umiConfig, final Statistics stats)
    {
        if(duplicateGroups == null)
            return;

        if(umiConfig.Enabled)
        {
            for(List<Fragment> fragments : duplicateGroups)
            {
                ++stats.DuplicateGroups;

                List<UmiGroup> umiGroups = buildUmiGroups(fragments, umiConfig);

                for(UmiGroup umiGroup : umiGroups)
                {
                    umiGroup.Fragments.forEach(x -> x.setUmiId(umiGroup.UmiId));

                    setDuplicateGroupProperties(umiGroup.Fragments);

                    ++stats.UmiGroups;
                    stats.addDuplicateGroup(umiGroup.Fragments);
                }
            }
        }
        else
        {
            for(List<Fragment> fragments : duplicateGroups)
            {
                setDuplicateGroupProperties(fragments);

                ++stats.DuplicateGroups;
                stats.addDuplicateGroup(fragments);
            }
        }
    }

    private static void setDuplicateGroupProperties(final List<Fragment> duplicateFragments)
    {
        int dupCount = duplicateFragments.size();
        duplicateFragments.forEach(x -> x.setDuplicateCount(dupCount));

        Fragment primary = findPrimaryFragment(duplicateFragments, false);
        primary.setStatus(PRIMARY);
    }

    public static Fragment findPrimaryFragment(final List<Fragment> fragments, boolean considerMarkedDups)
    {
        if(considerMarkedDups)
        {
            // take the primary (non-duplicate) group if there is (just) one already marked
            List<Fragment> nonDupGroups = fragments.stream().filter(x -> !hasDuplicates(x)).collect(Collectors.toList());

            if(nonDupGroups.size() == 1)
                return nonDupGroups.get(0);
        }

        // otherwise choose the group with the highest base quality
        Fragment maxFragment = null;
        double maxBaseQual = 0;

        for(Fragment fragment : fragments)
        {
            double avgBaseQual = calcBaseQualAverage(fragment);
            fragment.setAverageBaseQual(avgBaseQual);

            if(avgBaseQual > maxBaseQual)
            {
                maxBaseQual = avgBaseQual;
                maxFragment = fragment;
            }
        }

        return maxFragment;
    }

    private static boolean hasDuplicates(final Fragment fragment)
    {
        return fragment.reads().stream().anyMatch(x -> x.getDuplicateReadFlag());
    }

    public static double calcBaseQualAverage(final Fragment fragment)
    {
        int readBaseCount = 0;
        int readBaseQualTotal = 0;

        for(SAMRecord read : fragment.reads())
        {
            if(read.getSupplementaryAlignmentFlag())
                continue;

            for(int i = 0; i < read.getBaseQualities().length; ++i)
            {
                ++readBaseCount;
                readBaseQualTotal += read.getBaseQualities()[i];
            }
        }

        return readBaseCount > 0 ? readBaseQualTotal / (double)readBaseCount : 0;
    }

    public static List<UmiGroup> buildUmiGroups(final List<Fragment> fragments, final UmiConfig config)
    {
        Map<String,UmiGroup> groups = Maps.newHashMap();

        for(Fragment fragment : fragments)
        {
            String umiId = config.extractUmiId(fragment.id());

            UmiGroup group = groups.get(umiId);

            if(group == null)
            {
                groups.put(umiId, new UmiGroup(umiId, fragment));
            }
            else
            {
                group.Fragments.add(fragment);
            }
        }

        List<UmiGroup> orderedGroups = groups.values().stream().sorted(new SizeComparator()).collect(Collectors.toList());
        List<UmiGroup> finalGroups = Lists.newArrayList();

        int i = 0;
        while(i < orderedGroups.size() - 1)
        {
            UmiGroup first = orderedGroups.get(i);

            List<UmiGroup> cluster = Lists.newArrayList(first);

            int j = i + 1;
            while(j < orderedGroups.size())
            {
                UmiGroup second = orderedGroups.get(j);

                boolean merged = false;

                for(UmiGroup existing : cluster)
                {
                    if(existing.fragmentCount() >= second.fragmentCount() && !exceedsUmiIdDiff(existing.UmiId, second.UmiId))
                    {
                        merged = true;
                        break;
                    }
                }

                if(!merged)
                {
                    ++j;
                }
                else
                {
                    orderedGroups.remove(j);
                    cluster.add(second);

                    // restart the search since a newly added group may be close enough to a skipped one
                    j = i + 1;
                }
            }

            for(j = 1; j < cluster.size(); ++j)
            {
                first.Fragments.addAll(cluster.get(j).Fragments);
            }

            finalGroups.add(first);
            ++i;
        }

        return orderedGroups;
    }

    private static class SizeComparator implements Comparator<UmiGroup>
    {
        public int compare(final UmiGroup first, final UmiGroup second)
        {
            return first.Fragments.size() < second.Fragments.size() ? 1 : -1;
        }
    }


}
