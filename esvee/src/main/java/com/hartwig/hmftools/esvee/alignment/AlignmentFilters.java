package com.hartwig.hmftools.esvee.alignment;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_LOW_MOD_MQ_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ALIGNMENT_MIN_MOD_MAP_QUAL;
import static com.hartwig.hmftools.esvee.alignment.BreakendBuilder.segmentOrientation;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

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
            if(alignment.mapQual() >= ALIGNMENT_MIN_MOD_MAP_QUAL || alignment.hasAltAlignments())
            {
                candidateAlignments.add(alignment);
            }
            else
            {
                lowQualAlignments.add(alignment);
            }
        }

        // set modified map qual and then filtered low qual alignments
        Collections.sort(candidateAlignments, Comparator.comparingInt(x -> x.sequenceStart()));

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

        // first filter alignments with low modified map qual and no alt alignment info
        int validCount = (int)candidateAlignments.stream().filter(x -> x.exceedsMapQualThreshold()).count();

        if(candidateAlignments.size() == validCount)
        {
            validAlignments.addAll(candidateAlignments);
            return;
        }

        // for all the rest calculated an adjusted alignment score by subtracting overlap (inexact homology) and repeated bases from the score
        checkLocalVariants(candidateAlignments, validAlignments, lowQualAlignments);
    }

    private static void checkLocalVariants(
            final List<AlignData> candidateAlignments, final List<AlignData> validAlignments, final List<AlignData> lowQualAlignments)
    {
        if(candidateAlignments.size() < 2)
        {
            for(AlignData alignment : candidateAlignments)
            {
                if(alignment.exceedsMapQualThreshold())
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
            if(alignment.exceedsMapQualThreshold() == nextAlignment.exceedsMapQualThreshold())
                continue;

            findShortestLocalVariant(alignment, nextAlignment);
        }

        // now find linked short variants between adjacent alignments where both are below the required map qual threshold
        for(int i = 0; i < candidateAlignments.size() - 1; ++i)
        {
            AlignData alignment = candidateAlignments.get(i);

            AlignData nextAlignment = candidateAlignments.get(i + 1);

            // process any unlinked alignments
            if(alignment.exceedsMapQualThreshold() || nextAlignment.exceedsMapQualThreshold())
                continue;

            if(alignment.hasLinkedAltAlignment() && nextAlignment.hasLinkedAltAlignment())
                continue;

            findShortestLocalVariant(alignment, nextAlignment);
        }

        for(AlignData alignment : candidateAlignments)
        {
            if(alignment.exceedsMapQualThreshold() || alignment.hasLinkedAltAlignment() || alignment.useLowMapQualAlignment())
                validAlignments.add(alignment);
            else
                lowQualAlignments.add(alignment);
        }
    }

    private static boolean findShortestLocalVariant(final AlignData first, AlignData second)
    {
        // form alt alignments from the top alignment to make comparison easier
        List<AlternativeAlignment> firstAlignments, secondAlignments;
        AlternativeAlignment firstInitialAlignment, secondInitialAlignment;

        if(first.linkedAltAlignment() != null)
        {
            firstInitialAlignment = first.linkedAltAlignment();
            firstAlignments = List.of(firstInitialAlignment);
        }
        else
        {
            Orientation firstOrientation = segmentOrientation(first, true);
            int firstPosition = first.isForward() ? first.positionEnd() : first.positionStart();

            firstInitialAlignment = new AlternativeAlignment(
                    first.chromosome(), firstPosition, firstOrientation, "", first.mapQual());

            if(first.exceedsMapQualThreshold())
            {
                firstAlignments = List.of(firstInitialAlignment);
            }
            else
            {
                firstAlignments = Lists.newArrayList(firstInitialAlignment);
                firstAlignments.addAll(first.altAlignments());
            }
        }

        if(second.linkedAltAlignment() != null)
        {
            secondInitialAlignment = second.linkedAltAlignment();
            secondAlignments = List.of(secondInitialAlignment);
        }
        else
        {
            Orientation secondOrientation = segmentOrientation(second, false);
            int secondPosition = second.isForward() ? second.positionStart() : second.positionEnd();

            secondInitialAlignment = new AlternativeAlignment(
                    second.chromosome(), secondPosition, secondOrientation, "", second.mapQual());

            if(second.exceedsMapQualThreshold())
            {
                secondAlignments = List.of(secondInitialAlignment);
            }
            else
            {
                secondAlignments = Lists.newArrayList(secondInitialAlignment);
                secondAlignments.addAll(second.altAlignments());
            }
        }

        if(firstAlignments.size() < 2 && secondAlignments.size() < 2)
            return false;

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

        markAltAlignment(first, firstInitialAlignment, firstSelectedAlt);
        markAltAlignment(second, secondInitialAlignment, secondSelectedAlt);

        return firstSelectedAlt != null && secondSelectedAlt != null;
    }

    private static void markAltAlignment(
            final AlignData alignment, final AlternativeAlignment initialAlignment, final AlternativeAlignment altAlignment)
    {
        if(altAlignment == null)
            return;

        if(alignment.exceedsMapQualThreshold()) // only applicable for alignments failing the initial qual test
            return;

        if(initialAlignment == altAlignment) // allow the default alignment to be used
            alignment.markUseLowMapQualAlignment();
        else
            alignment.setLinkedAltAlignment(altAlignment);
    }

    @Deprecated
    @VisibleForTesting
    public static void filterAlignmentsOld(
            final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments, final List<AlignData> validAlignments,
            final List<AlignData> lowQualAlignments)
    {
        // exclude alignments with zero qual
        final List<AlignData> nonZeroAlignments = Lists.newArrayList();

        for(AlignData alignment : alignments)
        {
            if(alignment.mapQual() > 0)
                nonZeroAlignments.add(alignment);
            else
                lowQualAlignments.add(alignment);
        }

        if(nonZeroAlignments.isEmpty())
            return;

        // for all the rest calculated an adjusted alignment score by subtracting overlap (inexact homology) and repeated bases from the score
        String fullSequence = assemblyAlignment.fullSequence();

        nonZeroAlignments.forEach(x -> x.setFullSequenceData(fullSequence, assemblyAlignment.fullSequenceLength()));

        // set modified map qual and then filtered low qual alignments
        Collections.sort(nonZeroAlignments, Comparator.comparingInt(x -> x.sequenceStart()));

        for(int i = 0; i < nonZeroAlignments.size(); ++i)
        {
            AlignData alignment = nonZeroAlignments.get(i);

            int overlapStart = 0;
            int overlapEnd = 0;

            if(i > 0)
            {
                AlignData prevAlignment = nonZeroAlignments.get(i - 1);
                if(prevAlignment.sequenceEnd() >= alignment.sequenceStart())
                {
                    overlapStart = prevAlignment.sequenceEnd() - alignment.sequenceStart() + 1;
                }
            }

            if(i < nonZeroAlignments.size() - 1)
            {
                AlignData nextAlignment = nonZeroAlignments.get(i + 1);

                if(alignment.sequenceEnd() >= nextAlignment.sequenceStart())
                {
                    overlapEnd = alignment.sequenceEnd() - nextAlignment.sequenceStart() + 1;
                }
            }

            alignment.setAdjustedAlignment(fullSequence, overlapStart, overlapEnd);
        }

        for(AlignData alignment : nonZeroAlignments)
        {
            if(alignment.exceedsMapQualThreshold())
                validAlignments.add(alignment);
            else
                lowQualAlignments.add(alignment);
        }
    }
}
