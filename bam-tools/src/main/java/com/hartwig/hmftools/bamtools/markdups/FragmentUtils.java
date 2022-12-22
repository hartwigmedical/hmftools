package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import htsjdk.samtools.SAMRecord;

public class FragmentUtils
{
    public static int getUnclippedPosition(final SAMRecord read)
    {
        int position;

        if(orientation(read) == POS_ORIENT)
        {
            position = read.getAlignmentStart();
            if(read.getCigar().isLeftClipped())
                position -= read.getCigar().getFirstCigarElement().getLength();
        }
        else
        {
            position = read.getAlignmentEnd();
            if(read.getCigar().isRightClipped())
                position += read.getCigar().getLastCigarElement().getLength();
        }

        return position;
    }

    public static FragmentStatus calcFragmentStatus(final Fragment first, final Fragment second)
    {
        if(first.unpaired() != second.unpaired())
            return NONE;

        if(first.primaryReadsPresent() && second.primaryReadsPresent())
        {
            if(first.unpaired())
            {
                return first.initialPosition() == second.initialPosition() ? DUPLICATE : NONE;
            }
            else
            {
                return first.coordinates()[SE_START] == second.coordinates()[SE_START]
                        && first.coordinates()[SE_END] == second.coordinates()[SE_END] ? DUPLICATE : NONE;
            }
        }
        else
        {
            if(first.initialPosition() != second.initialPosition())
                return NONE;

            // mate start positions must be within close proximity
            SAMRecord firstRead = first.reads().get(0);
            SAMRecord secondRead = second.reads().get(0);

            if(!firstRead.getMateReferenceName().equals(secondRead.getMateReferenceName()))
                return NONE;

            return abs(firstRead.getMateAlignmentStart() - secondRead.getMateAlignmentStart()) < firstRead.getReadLength()
                    ? UNCLEAR : NONE;
        }
    }
    public static void classifyFragments(
            final PositionFragments positionFragments, final List<Fragment> resolvedFragments,
            final List<PositionFragments> incompletePositionFragments)
    {
        classifyFragments(positionFragments.Fragments, resolvedFragments, incompletePositionFragments);
    }

    public static void classifyFragments(
            final List<Fragment> positionFragments, final List<Fragment> resolvedFragments,
            final List<PositionFragments> incompletePositionFragments)
    {
        // take all the fragments at this initial fragment position and classify them as duplicates, non-duplicates or unclear

        // the list of fragments is copied and then reduced when duplicates or candidate duplicates are found
        List<Fragment> allFragments = Lists.newArrayList(positionFragments);

        if(allFragments.size() == 1)
        {
            Fragment fragment = allFragments.get(0);
            fragment.setStatus(NONE);
            resolvedFragments.add(fragment);
            return;
        }

        int i = 0;
        while(i < allFragments.size())
        {
            Fragment fragment1 = allFragments.get(i);

            if(i == allFragments.size() - 1)
            {
                fragment1.setStatus(NONE);
                resolvedFragments.add(fragment1);
                break;
            }

            PositionFragments incompleteFragments = null;
            List<Fragment> duplicateFragments = null;

            int j = i + 1;

            while(j < allFragments.size())
            {
                Fragment fragment2 = allFragments.get(j);

                FragmentStatus status = calcFragmentStatus(fragment1, fragment2);

                if(status == NONE)
                {
                    ++j;
                    continue;
                }

                if(fragment1.status() != UNSET && status != NONE && fragment1.status() != status)
                {
                    BM_LOGGER.warn("fragment({}) has alt status({}) with other({})", fragment1, status, fragment2);
                }
                else if(status == DUPLICATE)
                {
                    fragment1.setStatus(status);
                    fragment2.setStatus(status);

                    if(duplicateFragments == null)
                        duplicateFragments = Lists.newArrayList(fragment1);

                    duplicateFragments.add(fragment2);
                    allFragments.remove(j);
                }
                else if(status == UNCLEAR)
                {
                    fragment1.setStatus(status);
                    fragment2.setStatus(status);

                    if(incompleteFragments == null)
                        incompleteFragments = new PositionFragments(fragment1.initialPosition(), fragment1);

                    incompleteFragments.Fragments.add(fragment2);
                    allFragments.remove(j);
                }
            }

            if(fragment1.status().isDuplicate())
            {
                resolvedFragments.addAll(duplicateFragments);

                Fragment primary = findPrimaryFragment(duplicateFragments, true);
                primary.setStatus(PRIMARY);
            }
            else if(incompleteFragments != null)
            {
                incompletePositionFragments.add(incompleteFragments);
            }
            else
            {
                fragment1.setStatus(NONE);
                resolvedFragments.add(fragment1);
            }

            ++i;
        }
    }

    private static boolean hasDuplicates(final Fragment fragment)
    {
        return fragment.reads().stream().anyMatch(x -> x.getDuplicateReadFlag());
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
        int maxBaseQual = 0;

        for(Fragment fragment : fragments)
        {
            int groupBaseQual = calcBaseQualTotal(fragment);

            if(groupBaseQual > maxBaseQual)
            {
                maxBaseQual = groupBaseQual;
                maxFragment = fragment;
            }
        }

        return maxFragment;
    }

    public static int calcBaseQualTotal(final Fragment fragment)
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

        return readBaseCount > 0 ? (int)round(readBaseQualTotal / (double)readBaseCount) : 0;
    }

    private static final String CHR_PARTITION_DELIM = "_";

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosome + CHR_PARTITION_DELIM + partition;
    }

    @Deprecated
    public static void reconcileFragments(
            final Map<String,Fragment> supplementaries, final Map<String,Fragment> resolvedFragments,
            final List<PositionFragments> incompletePositionFragments)
    {
        // first add any supplementaries or incomplete fragments to resolved fragments

        // link up any fragments by read ID and look for complete fragments
        Set<PositionFragments> modifiedPositionFragments = Sets.newHashSet();

        Map<String,Fragment> incompleteFragments = Maps.newHashMap();
        incompletePositionFragments.forEach(x -> x.Fragments.forEach(y -> incompleteFragments.put(y.id(), y)));

        for(Fragment fragment : resolvedFragments.values())
        {
            Fragment supp = supplementaries.get(fragment.id());

            if(supp != null)
            {
                supp.reads().forEach(x -> fragment.addRead(x));
                supplementaries.remove(supp.id());
            }

            Fragment incompleteFrag = incompleteFragments.get(fragment.id());

            if(incompleteFrag != null)
            {
                incompleteFrag.reads().forEach(x -> fragment.addRead(x));
                incompleteFragments.remove(supp.id());

                PositionFragments positionFragments = incompletePositionFragments.stream()
                        .filter(x -> x.Position == incompleteFrag.initialPosition()).findFirst().orElse(null);

                if(positionFragments != null)
                {
                    positionFragments.Fragments.remove(incompleteFrag);
                    modifiedPositionFragments.add(positionFragments);
                }
            }
        }

        for(PositionFragments positionFragments : modifiedPositionFragments)
        {
            if(positionFragments.Fragments.size() < 2)
            {
                for(Fragment fragment : positionFragments.Fragments)
                {
                    fragment.setStatus(NONE);
                    resolvedFragments.put(fragment.id(), fragment);
                }

                incompletePositionFragments.remove(positionFragments);
            }
        }
    }

    /*

                if(mConfig.RunChecks)
            {
                // log discrepancies
                ReadGroup calcPrimaryGroup = duplicateGroup.findPrimaryGroup(false);

                boolean logDiscrepancy = false;

                if(primaryGroup != calcPrimaryGroup && calcBaseQualTotal(primaryGroup) != calcBaseQualTotal(calcPrimaryGroup))
                    logDiscrepancy = true;
                else if(duplicateGroup.readGroups().stream().anyMatch(x -> hasDuplicates(x) == (x == primaryGroup)))
                    logDiscrepancy = true;

                if(logDiscrepancy)
                {
                    for(ReadGroup readGroup : duplicateGroup.readGroups())
                    {
                        BM_LOGGER.trace("readGroup({}) hasDups({}) isPrimary({}) baseQualTotal({})",
                                readGroup.toString(), hasDuplicates(readGroup), readGroup == primaryGroup, calcBaseQualTotal(readGroup));
                    }
                }
            }
     */

}