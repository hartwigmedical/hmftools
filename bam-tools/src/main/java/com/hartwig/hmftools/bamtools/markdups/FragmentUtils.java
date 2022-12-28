package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNCLEAR;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

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


    // initial call needs to prioritise finding resolved fragments (dups or NONE), and form unclear groups from amongst the rest


    // second call need only look for resolved fragments

    public static void classifyFragments(
            final List<Fragment> fragments, final List<Fragment> resolvedFragments,
            @Nullable final List<CandidateDuplicates> candidateDuplicatesList)
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

        // at most 1 position with unclear fragments will be created since having more can lead to overlapping groups

        int fragmentCount = fragments.size();
        Set<Fragment> possibleDuplicates = Sets.newHashSet();

        int i = 0;
        while(i < fragments.size())
        {
            Fragment fragment1 = fragments.get(i);

            if(i == fragments.size() - 1)
            {
                if(!possibleDuplicates.contains(fragment1))
                {
                    fragment1.setStatus(NONE);
                    resolvedFragments.add(fragment1);
                    fragments.remove(i);
                }
                break;
            }

            List<Fragment> duplicateFragments = null;

            int j = i + 1;
            while(j < fragments.size())
            {
                Fragment fragment2 = fragments.get(j);

                FragmentStatus status = calcFragmentStatus(fragment1, fragment2);

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

                if(fragment1.status() != DUPLICATE && status == UNCLEAR)
                {
                    // the pair is a candidate for duplicates but without their mates it's unclear whether they will be
                    possibleDuplicates.add(fragment1);
                    possibleDuplicates.add(fragment2);
                }

                ++j;
            }

            if(fragment1.status().isDuplicate())
            {
                resolvedFragments.addAll(duplicateFragments);
                fragments.remove(i);

                Fragment primary = findPrimaryFragment(duplicateFragments, true);
                primary.setStatus(PRIMARY);
            }
            else if(possibleDuplicates.contains(fragment1))
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

        List<Fragment> unclearFragments = possibleDuplicates.stream().filter(y -> !resolvedFragments.contains(y)).collect(Collectors.toList());
        unclearFragments.forEach(x -> x.setStatus(UNCLEAR));

        if(candidateDuplicatesList != null && !unclearFragments.isEmpty())
        {
            candidateDuplicatesList.add(new CandidateDuplicates(unclearFragments.get(0).initialPosition(), unclearFragments));
        }

        if(unclearFragments.size() + resolvedFragments.size() != fragmentCount)
        {
            BM_LOGGER.error("failed to classify all fragments: original({}) resolved({}) unclear({})",
                    fragmentCount, resolvedFragments.size(), unclearFragments.size());
        }
    }

    public static void checkFragmentClassification(
            final List<Fragment> resolvedFragments, final List<CandidateDuplicates> incompletePositionFragments)
    {
        if(resolvedFragments.stream().anyMatch(x -> x.status() == UNSET))
        {
            BM_LOGGER.error("failed to classify all resolved fragments");
        }

        if(incompletePositionFragments.stream().anyMatch(x -> x.Fragments.stream().anyMatch(y -> y.status() == UNSET)))
        {
            BM_LOGGER.error("failed to classify all incomplete fragments");
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
        return chromosomeIndicator(chromosome) + partition;
    }

    public static String chromosomeIndicator(final String chromosome)
    {
        return chromosome + CHR_PARTITION_DELIM;
    }

    public static String readToString(final SAMRecord read)
    {
        return format("id(%s) coords(coords(%s:%d-%d) cigar(%s) mate(%s:%d) flags(%d)",
                read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                read.getCigarString(), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getFlags());
    }

    public static boolean readInSpecifiedRegions(
            final SAMRecord read, final List<ChrBaseRegion> regions, final List<String> chromosomes)
    {
        if(!chromosomes.isEmpty())
        {
            if(chromosomes.stream().noneMatch(x -> x.equals(read.getContig())))
                return false;

            // any mates or supplementaries must also be within the regions specified
            if(chromosomes.stream().noneMatch(x -> x.equals(read.getMateReferenceName())))
                return false;
        }

        if(!regions.isEmpty())
        {
            if(regions.stream().noneMatch(x -> x.containsPosition(read.getContig(), read.getAlignmentStart())))
                return false;

            // any mates or supplementaries must also be within the regions specified
            if(regions.stream().noneMatch(x -> x.containsPosition(read.getMateReferenceName(), read.getMateAlignmentStart())))
                return false;
        }

        // ignore checking supplementaries since a) they aren't marked as duplicates by other tools and b) they shouldn't be a reason
        // to ignore a primary read since that then impacts duplicate classification
        /*
        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            SupplementaryReadData suppData = SupplementaryReadData.from(read);

            if(suppData != null)
            {
                if(!regions.isEmpty() && regions.stream().noneMatch(x -> x.containsPosition(suppData.Chromosome, suppData.Position)))
                    return false;

                if(!chromosomes.isEmpty() && chromosomes.stream().noneMatch(x -> x.equals(suppData.Chromosome)))
                    return false;
            }
        }
        */

        return true;
    }
}