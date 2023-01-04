package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentCoordinates.NO_COORDS;
import static com.hartwig.hmftools.bamtools.markdups.FragmentCoordinates.formCoordinate;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

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

    public static FragmentCoordinates getFragmentCoordinates(final List<SAMRecord> reads)
    {
        SAMRecord firstRead = null;
        SAMRecord mateRead = null;

        for(SAMRecord read : reads)
        {
            if(read.getSupplementaryAlignmentFlag())
                continue;

            if(firstRead == null)
            {
                firstRead = read;
            }
            else
            {
                mateRead = read;
                break;
            }
        }

        boolean readForwardStrand = orientation(firstRead) == POS_ORIENT;

        int readCoordinate = firstRead.getCigar() != null ?
                getUnclippedPosition(firstRead) : getUnclippedPosition(firstRead.getAlignmentStart(), firstRead.getCigarString(), readForwardStrand);

        int readStrandPosition = readForwardStrand ? readCoordinate : -readCoordinate;
        String readCoordStr = formCoordinate(firstRead.getReferenceName(), readCoordinate, readForwardStrand);

        if(!firstRead.getReadPairedFlag())
            return new FragmentCoordinates(readCoordStr, readStrandPosition);

        if(mateRead == null)
            return new FragmentCoordinates(readCoordStr, readStrandPosition, true);

        boolean mateForwardStrand = orientation(mateRead) == POS_ORIENT;

        int mateCoordinate = mateRead.getCigar() != null ?
                getUnclippedPosition(mateRead) : getUnclippedPosition(mateRead.getAlignmentStart(), mateRead.getCigarString(), mateForwardStrand);

        int mateStrandPosition = mateForwardStrand ? mateCoordinate : -mateCoordinate;
        String mateCoordStr = formCoordinate(mateRead.getReferenceName(), mateCoordinate, mateForwardStrand);

        boolean readLowerPos;
        if(firstRead.getReferenceIndex() == firstRead.getMateReferenceIndex())
        {
            readLowerPos = readCoordinate <= mateCoordinate;
        }
        else
        {
            readLowerPos = firstRead.getReferenceIndex() < firstRead.getMateReferenceIndex();
        }

        return readLowerPos ?
                new FragmentCoordinates(readCoordStr + "_" + mateCoordStr, readStrandPosition)
                : new FragmentCoordinates(mateCoordStr + "_" + readCoordStr, mateStrandPosition);
    }

    // public static final int PROXIMATE_MATE_INSERT_DIFF = 100;

    public static FragmentStatus calcFragmentStatus(final Fragment first, final Fragment second)
    {
        if(first.unpaired() != second.unpaired())
            return NONE;

        if(first.primaryReadsPresent() && second.primaryReadsPresent())
            return first.coordinates().Key.equals(second.coordinates().Key) ? DUPLICATE : NONE;

        if(first.initialPosition() != second.initialPosition())
            return NONE;

        // mate start positions must be within close proximity
        SAMRecord firstRead = first.reads().get(0);
        SAMRecord secondRead = second.reads().get(0);

        if(!firstRead.getMateReferenceName().equals(secondRead.getMateReferenceName()))
            return NONE;

        if(firstRead.getMateNegativeStrandFlag() != secondRead.getMateNegativeStrandFlag())
            return NONE;

        return abs(firstRead.getMateAlignmentStart() - secondRead.getMateAlignmentStart()) < firstRead.getReadLength()
                ? CANDIDATE : NONE;
    }

    private static final int HIGH_DEPTH_THRESHOLD = 10000;

    public static void classifyFragments(
            final List<Fragment> fragments, final List<Fragment> resolvedFragments,
            final List<CandidateDuplicates> candidateDuplicatesList)
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
        boolean applyHighDepthLogic = fragmentCount > HIGH_DEPTH_THRESHOLD;

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

                Fragment primary = findPrimaryFragment(duplicateFragments, true);
                primary.setStatus(PRIMARY);

                // apply UMI logic and create a consensus read here

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
        return format("id(%s) coords(%s:%d-%d) cigar(%s) mate(%s:%d) flags(%d)",
                read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                read.getCigarString(), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getFlags());
    }

    public static boolean readInSpecifiedRegions(
            final SAMRecord read, final List<ChrBaseRegion> regions, final List<String> chromosomes, boolean checkSupplementaries)
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

        // by default ignore checking supplementaries since a) they aren't marked as duplicates by other tools and b) they shouldn't be
        // a reason to ignore a primary read since that then impacts duplicate classification
        if(checkSupplementaries)
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

        return true;
    }

    // methods which are reliant on having mate CIGAR:
    public static FragmentCoordinates getFragmentCoordinates(final SAMRecord read, boolean orderCoordinates)
    {
        boolean readForwardStrand = orientation(read) == POS_ORIENT;

        int readCoordinate = read.getCigar() != null ?
                getUnclippedPosition(read) : getUnclippedPosition(read.getAlignmentStart(), read.getCigarString(), readForwardStrand);

        int readStrandPosition = readForwardStrand ? readCoordinate : -readCoordinate;
        String readCoordStr = formCoordinate(read.getReferenceName(), readCoordinate, readForwardStrand);

        if(!read.getReadPairedFlag())
            return new FragmentCoordinates(readCoordStr, readStrandPosition);

        if(!read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
            return NO_COORDS;

        String mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        boolean mateForwardStrand = !read.getMateNegativeStrandFlag();
        int mateCoordinate = getUnclippedPosition(read.getMateAlignmentStart(), mateCigar, mateForwardStrand);
        int mateStrandPosition = mateForwardStrand ? mateCoordinate : -mateCoordinate;

        String mateCoordStr = formCoordinate(read.getMateReferenceName(), mateCoordinate, mateForwardStrand);

        if(!orderCoordinates)
            return new FragmentCoordinates(readCoordStr + "_" + mateCoordStr, readStrandPosition);

        boolean readLowerPos;
        if(read.getReferenceIndex() == read.getMateReferenceIndex())
        {
            readLowerPos = readCoordinate <= mateCoordinate;
        }
        else
        {
            readLowerPos = read.getReferenceIndex() < read.getMateReferenceIndex();
        }

        return readLowerPos ?
                new FragmentCoordinates(readCoordStr + "_" + mateCoordStr, readStrandPosition)
                : new FragmentCoordinates(mateCoordStr + "_" + readCoordStr, mateStrandPosition);
    }

    public static int getUnclippedPosition(final int readStart, final String cigarStr, final boolean forwardStrand)
    {
        int currentPosition = readStart;
        int elementLength = 0;

        for(int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            boolean isAddItem = (c == 'D' || c == 'M' || c == 'S' || c == 'N');

            if(isAddItem)
            {
                if(forwardStrand)
                {
                    // back out the left clip if present
                    return c == 'S' ? readStart - elementLength : readStart;
                }

                if(c == 'S' && readStart == currentPosition)
                {
                    // ignore left-clip when getting reverse strand position
                }
                else
                {
                    currentPosition += elementLength;
                }

                elementLength = 0;
                continue;
            }

            int digit = c - '0';
            if (digit >= 0 && digit <= 9)
            {
                elementLength = elementLength * 10 + digit;
            }
            else
            {
                elementLength = 0;
            }
        }

        // always pointing to the start of the next element, so need to move back a base
        return currentPosition - 1;
    }
}