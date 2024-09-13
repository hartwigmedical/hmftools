package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.bam.CigarUtils.calcCigarAlignedLength;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_LOW_MOD_MQ_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_ADJUST_ALIGN_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_MOD_MAP_QUAL;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_MOD_MAP_QUAL_NO_XA;
import static com.hartwig.hmftools.esvee.alignment.BreakendBuilder.segmentOrientation;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;

public final class AlignmentFilters
{
    @VisibleForTesting
    public static void filterAlignments(
            final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments, final List<AlignData> validAlignments,
            final List<AlignData> lowQualAlignments)
    {
        // set sequence attributes and then modified map qual
        String fullSequence = assemblyAlignment.fullSequence();

        alignments.forEach(x -> x.setFullSequenceData(fullSequence, assemblyAlignment.fullSequenceLength()));

        // remove low-qual alignments with no alternatives
        List<AlignData> candidateAlignments = Lists.newArrayList();

        for(AlignData alignment : alignments)
        {
            if(exceedsInitialMapQualThreshold(alignment) || alignment.hasAltAlignments())
            {
                candidateAlignments.add(alignment);
            }
            else
            {
                lowQualAlignments.add(alignment);
            }
        }

        // sort alignments but leave any matching sequence starts as they are
        Collections.sort(candidateAlignments, new AlignmentOrderComparator());

        // set modified map qual and then filtered low qual alignments
        for(int i = 0; i < candidateAlignments.size(); ++i)
        {
            AlignData alignment = candidateAlignments.get(i);

            int overlapStart = 0;
            int overlapEnd = 0;

            if(i > 0)
            {
                AlignData prevAlignment = candidateAlignments.get(i - 1);
                if(prevAlignment.sequenceEnd() >= alignment.sequenceStart())
                {
                    overlapStart = prevAlignment.sequenceEnd() - alignment.sequenceStart() + 1;
                }
            }

            if(i < candidateAlignments.size() - 1)
            {
                AlignData nextAlignment = candidateAlignments.get(i + 1);

                if(alignment.sequenceEnd() >= nextAlignment.sequenceStart())
                {
                    overlapEnd = alignment.sequenceEnd() - nextAlignment.sequenceStart() + 1;
                }
            }

            alignment.setAdjustedAlignment(fullSequence, overlapStart, overlapEnd);
        }

        // check adjusted alignment lengths for inner alignments  only
        int index = 1;
        while(index < candidateAlignments.size() - 1)
        {
            AlignData alignment = candidateAlignments.get(index);

            // inner alignments must exceed min aligned length
            if(alignment.adjustedAlignment() < ALIGNMENT_MIN_ADJUST_ALIGN_LENGTH)
            {
                candidateAlignments.remove(index);
                lowQualAlignments.add(alignment);
            }
            else
            {
                ++index;
            }
        }

        // check modified map quals
        index = 0;
        while(index < candidateAlignments.size())
        {
            AlignData alignment = candidateAlignments.get(index);

            if(!exceedsNoAltsModMapQualThreshold(alignment) && !alignment.hasAltAlignments())
            {
                candidateAlignments.remove(index);
                lowQualAlignments.add(alignment);
            }
            else
            {
                ++index;
            }
        }
        int validCount = (int)candidateAlignments.stream().filter(x -> exceedsModMapQualThreshold(x)).count();

        if(candidateAlignments.size() == validCount)
        {
            validAlignments.addAll(candidateAlignments);
            return;
        }

        // for all the rest calculated an adjusted alignment score by subtracting overlap (inexact homology) and repeated bases from the score
        checkLocalVariants(candidateAlignments, validAlignments, lowQualAlignments);
    }

    private static boolean exceedsInitialMapQualThreshold(final AlignData alignment)
    {
        return exceedsMapQualThreshold(alignment, false, ALIGNMENT_MIN_MOD_MAP_QUAL_NO_XA);
    }

    private static boolean exceedsModMapQualThreshold(final AlignData alignment)
    {
        return exceedsMapQualThreshold(alignment, true, ALIGNMENT_MIN_MOD_MAP_QUAL);
    }

    private static boolean exceedsNoAltsModMapQualThreshold(final AlignData alignment)
    {
        return exceedsMapQualThreshold(alignment, true, ALIGNMENT_MIN_MOD_MAP_QUAL_NO_XA);
    }

    private static boolean exceedsMapQualThreshold(final AlignData alignment, boolean useModified, int threshold)
    {
        double mapQual = useModified ? alignment.modifiedMapQual() : alignment.mapQual();
        return mapQual >= threshold;
    }

    private static class AlignmentOrderComparator implements Comparator<AlignData>
    {
        @Override
        public int compare(final AlignData first, final AlignData second)
        {
            if(first.sequenceStart() == second.sequenceStart())
                return -1;

            return first.sequenceStart() < second.sequenceStart() ? -1 : 1;
        }
    }

    private static void checkLocalVariants(
            final List<AlignData> candidateAlignments, final List<AlignData> validAlignments, final List<AlignData> lowQualAlignments)
    {
        // first exit if there aren't enough alignments to consider making links using those of lower map qual
        if(candidateAlignments.size() < 2)
        {
            for(AlignData alignment : candidateAlignments)
            {
                if(exceedsModMapQualThreshold(alignment))
                    validAlignments.add(alignment);
                else
                    lowQualAlignments.add(alignment);
            }

            return;
        }

        // first mark candidates with sufficient mod map qual and find linked short variants with those
        for(int i = 0; i < candidateAlignments.size() - 1; ++i)
        {
            AlignData alignment = candidateAlignments.get(i);

            AlignData nextAlignment = candidateAlignments.get(i + 1);

            // process if one valid and the other not
            if(exceedsModMapQualThreshold(alignment) == exceedsModMapQualThreshold(nextAlignment))
                continue;

            findShortestLocalVariant(alignment, nextAlignment);
        }

        // now find linked short variants between adjacent alignments where both are below the required map qual threshold
        for(int i = 0; i < candidateAlignments.size() - 1; ++i)
        {
            AlignData alignment = candidateAlignments.get(i);

            AlignData nextAlignment = candidateAlignments.get(i + 1);

            // process any unlinked alignments
            if(exceedsModMapQualThreshold(alignment) || exceedsModMapQualThreshold(nextAlignment))
                continue;

            if(alignment.hasLowMapQualAlignment() && nextAlignment.hasLowMapQualAlignment())
                continue;

            findShortestLocalVariant(alignment, nextAlignment);
        }

        for(AlignData alignment : candidateAlignments)
        {
            if(exceedsModMapQualThreshold(alignment) || alignment.hasLowMapQualAlignment())
                validAlignments.add(alignment);
            else
                lowQualAlignments.add(alignment);
        }
    }

    private static List<AlternativeAlignment> collectAlignmentCandidates(final AlignData alignment, boolean linksAtEnd)
    {
        List<AlternativeAlignment> alignments;

        if(alignment.selectedAltAlignment() != null)
        {
            AlternativeAlignment initialAlignment = alignment.selectedAltAlignment();
            alignments = List.of(initialAlignment);
        }
        else
        {
            Orientation orientation = segmentOrientation(alignment, linksAtEnd);

            // links at end, meaning this segment's end is linked to the start of another, so use its end position if not reversed
            boolean useEndPosition = alignment.isForward() == linksAtEnd;
            int position = useEndPosition ? alignment.positionEnd() : alignment.positionStart();

            AlternativeAlignment initialAlignment = new AlternativeAlignment(
                    alignment.chromosome(), position, orientation, "", alignment.mapQual());

            if(exceedsModMapQualThreshold(alignment))
            {
                alignments = List.of(initialAlignment);
            }
            else
            {
                alignments = Lists.newArrayList(initialAlignment);

                for(AlternativeAlignment altAlignment : alignment.rawAltAlignments())
                {
                    // apply orientation info
                    Orientation altOrientation = segmentOrientation(altAlignment.Orient, linksAtEnd);

                    useEndPosition = altAlignment.Orient.isForward() == linksAtEnd;
                    position = altAlignment.Position;

                    if(useEndPosition)
                        position += calcCigarAlignedLength(altAlignment.Cigar) - 1;

                    alignments.add(new AlternativeAlignment(
                            altAlignment.Chromosome,  position, altOrientation, altAlignment.Cigar, altAlignment.MapQual));
                }
            }
        }

        return alignments;
    }

    private static boolean findShortestLocalVariant(final AlignData first, AlignData second)
    {
        // form alt alignments from the top alignment to make comparison easier
        List<AlternativeAlignment> firstAlignments = collectAlignmentCandidates(first, true);
        List<AlternativeAlignment> secondAlignments = collectAlignmentCandidates(second, false);

        if(firstAlignments.isEmpty() || secondAlignments.isEmpty())
            return false;

        AlternativeAlignment firstInitialAlignment = firstAlignments.get(0);
        AlternativeAlignment secondInitialAlignment = secondAlignments.get(0);

        int shortestLength = -1;
        AlternativeAlignment firstSelectedAlt = null;
        AlternativeAlignment secondSelectedAlt = null;

        for(AlternativeAlignment firstAlt : firstAlignments)
        {
            for(AlternativeAlignment secondAlt : secondAlignments)
            {
                if(!firstAlt.Chromosome.equals(secondAlt.Chromosome))
                    continue;

                int svLength = abs(firstAlt.Position - secondAlt.Position);

                if(svLength > ALIGNMENT_LOW_MOD_MQ_VARIANT_LENGTH)
                    continue;

                if(shortestLength < 0 || svLength < shortestLength)
                {
                    firstSelectedAlt = firstAlt;
                    secondSelectedAlt = secondAlt;
                    shortestLength = svLength;
                }
            }
        }

        if(firstSelectedAlt == null)
            firstSelectedAlt = firstInitialAlignment;

        if(secondSelectedAlt == null)
            secondSelectedAlt = secondInitialAlignment;

        boolean hasShortSvLink = shortestLength >= 0;
        markAltAlignment(first, firstInitialAlignment, firstSelectedAlt, firstAlignments, hasShortSvLink);
        markAltAlignment(second, secondInitialAlignment, secondSelectedAlt, secondAlignments, hasShortSvLink);

        return firstSelectedAlt != null && secondSelectedAlt != null;
    }

    private static void markAltAlignment(
            final AlignData alignment, final AlternativeAlignment initialAlignment, final AlternativeAlignment selectedAlignment,
            final List<AlternativeAlignment> allAlignments, boolean hasShortSvLink)
    {
        if(selectedAlignment == null)
            return;

        if(exceedsModMapQualThreshold(alignment)) // only applicable for alignments failing the initial qual test
            return;

        // by marking either the selected alt alignment or registered alternatives, the low map-qual alignment can now be used as a breakend
        List<AlternativeAlignment> unselectedAltAlignments = allAlignments.stream()
                .filter(x -> x != selectedAlignment).collect(Collectors.toList());

        alignment.setSelectedAltAlignments(
                selectedAlignment != initialAlignment ? selectedAlignment : null, unselectedAltAlignments, hasShortSvLink);
    }
}
