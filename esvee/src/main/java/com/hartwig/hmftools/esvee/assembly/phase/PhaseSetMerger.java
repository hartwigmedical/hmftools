package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createMinBaseQuals;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.SequenceCompare;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

public final class PhaseSetMerger
{
    public static void mergePhaseSets(final List<PhaseSet> phaseSets)
    {
        if(phaseSets.size() == 1)
            return;

        for(int i = 0; i < phaseSets.size() - 1; ++i)
        {
            PhaseSet first = phaseSets.get(i);

            if(!first.assemblyAlignment().isValid())
                continue;

            for(int j = i + 1; j < phaseSets.size(); ++j)
            {
                PhaseSet second = phaseSets.get(j);

                if(!second.assemblyAlignment().isValid())
                    continue;

                tryMergePhaseSets(first, second);
            }
        }
    }

    public static void tryMergePhaseSets(final PhaseSet firstPhaseSet, final PhaseSet secondPhaseSet)
    {
        AssemblyAlignment firstAssembly = firstPhaseSet.assemblyAlignment();
        AssemblyAlignment secondAssembly = secondPhaseSet.assemblyAlignment();

        boolean firstIsPrimary = firstAssembly.fullSequenceLength() >= secondAssembly.fullSequenceLength();

        AssemblyAlignment primaryAssembly, otherAssembly;

        // PhaseSet primaryPhaseSet = selectPrimaryPhaseSet(firstPhaseSet, secondPhaseSet);

        PhaseSet primaryPhaseSet, otherPhaseSet;

        if(firstIsPrimary)
        {
            primaryAssembly = firstAssembly;
            otherAssembly = secondAssembly;
            primaryPhaseSet = firstPhaseSet;
            otherPhaseSet = secondPhaseSet;
        }
        else
        {
            primaryAssembly = secondAssembly;
            otherAssembly = firstAssembly;
            primaryPhaseSet = secondPhaseSet;
            otherPhaseSet = firstPhaseSet;
        }

        String primarySequence = primaryAssembly.fullSequence();
        String otherFullSequence = otherAssembly.fullSequence();

        for(int i = 0; i <= 1; ++i)
        {
            Orientation otherOrientation = (i == 0) ? REVERSE : FORWARD;

            String otherSequence = otherOrientation == FORWARD ?
                    otherFullSequence : Nucleotides.reverseComplementBases(otherFullSequence);

            int[] matchIndices = findMergeIndices(primaryAssembly.fullSequence(), otherSequence);

            if(matchIndices != null)
            {
                mergeAssemblyAlignments(primaryAssembly, otherAssembly, primarySequence, otherSequence, otherOrientation, matchIndices);

                primaryPhaseSet.mergePhaseSet(otherPhaseSet);
                otherPhaseSet.setMergedPhaseSetId(primaryPhaseSet.id());

                return;
            }
        }
    }

    public static final int MATCH_SUBSEQUENCE_LENGTH = 50;

    private static int[] findMergeIndices(final String firstSequence, final String secondSequence)
    {
        // if a match is found, returns the first and seconds' match start indices
        // one or the other will have a value of zero, meaning the start of its sequence matches within the other

        // first attempt a simple (optimistic) match
        int matchIndex = firstSequence.indexOf(secondSequence);

        if(matchIndex >= 0)
            return new int[] { matchIndex, 0 };

        // try subsequence matches
        int firstLength = firstSequence.length();
        int secondLength = secondSequence.length();

        List<RepeatInfo> firstRepeats = null;
        List<RepeatInfo> secondRepeats = null;

        byte[] firstBases = null;
        byte[] firstBaseQuals = null;
        byte[] secondBases = null;
        byte[] secondBaseQuals = null;

        int matchSeqStartIndex = 0;

        int subSeqIterations = (int)floor(firstLength / MATCH_SUBSEQUENCE_LENGTH);

        for(int i = 0; i < subSeqIterations; ++i) // being the total junction sequence length (ie 100) divided by the subsequence length
        {
            matchSeqStartIndex = i * MATCH_SUBSEQUENCE_LENGTH;
            int matchSeqEndIndex = matchSeqStartIndex + MATCH_SUBSEQUENCE_LENGTH;

            if(matchSeqEndIndex >= firstLength + 1)
                break;

            String firstSubSequence = firstSequence.substring(matchSeqStartIndex, matchSeqEndIndex);

            int secondSubSeqIndex = secondSequence.indexOf(firstSubSequence);

            if(secondSubSeqIndex < 0)
                continue;

            // from this subsquence, expand to cover the maximal overlap of the two sequences
            int minStartDistance = min(matchSeqStartIndex, secondSubSeqIndex);

            int firstMatchIndexStart = matchSeqStartIndex - minStartDistance;
            int secondMatchIndexStart = secondSubSeqIndex - minStartDistance;

            int overlapLength = min(firstLength - firstMatchIndexStart - 1, secondLength - secondMatchIndexStart - 1);

            int firstMatchIndexEnd = overlapLength + 1;
            int secondMatchIndexEnd = overlapLength + 1;

            if(firstBases == null && secondBases == null)
            {
                firstBases = firstSequence.getBytes();
                firstBaseQuals = createMinBaseQuals(firstLength);
                firstRepeats = RepeatInfo.findRepeats(firstBases);

                secondBases = secondSequence.getBytes();
                secondBaseQuals = createMinBaseQuals(secondLength);
                secondRepeats = RepeatInfo.findRepeats(secondBases);
            }

            int permittedMismatches = mismatchesPerComparisonLength(overlapLength);

            int mismatchCount = SequenceCompare.compareSequences(
                    firstBases, firstBaseQuals, firstMatchIndexStart, firstMatchIndexEnd, firstRepeats,
                    secondBases, secondBaseQuals, secondMatchIndexStart, secondMatchIndexEnd, secondRepeats, permittedMismatches);

            if(mismatchCount <= permittedMismatches)
            {
                return new int[] { firstMatchIndexStart, secondMatchIndexStart };
            }
        }

        return null;
    }

    public static void mergeAssemblyAlignments(
            final AssemblyAlignment primaryAssembly, final AssemblyAlignment otherAssembly,
            final String primarySequence, final String otherSequence, final Orientation otherOrientation, int[] matchIndices)
    {
        // first merge the sequences
        int primaryMatchStart = matchIndices[0];
        int otherMatchStart = matchIndices[1];

        int primaryLength = primarySequence.length();
        int otherLength = otherAssembly.fullSequenceLength();

        String newSequence = null;

        int primaryOffsetAdjust = 0;
        int otherOffsetAdjust = 0;
        Map<Integer,String> newSequenceOverlaps = Maps.newHashMap();

        if(otherMatchStart > primaryMatchStart)
        {
            // the non-primary assemblie starts earlier than the primrary so extend at the start
            primaryOffsetAdjust = otherMatchStart - primaryMatchStart;
            newSequence = otherSequence.substring(0, primaryOffsetAdjust);
            newSequence += primarySequence;

            int minOtherSeqIndex = primaryOffsetAdjust;
            otherAssembly.sequenceOverlaps().entrySet().stream()
                    .filter(x -> x.getKey() < minOtherSeqIndex).forEach(x -> newSequenceOverlaps.put(x.getKey(), x.getValue()));
        }
        else
        {
            // the primary starts earlier
            otherOffsetAdjust = primaryMatchStart - otherMatchStart;
        }

        if(primaryMatchStart + otherLength > primaryLength)
        {
            if(newSequence == null)
                newSequence = primarySequence;

            int extensionLength = primaryMatchStart + otherLength - primaryLength;
            int otherSeqIndex = otherLength - extensionLength;

            newSequence += otherSequence.substring(otherSeqIndex);

            otherAssembly.sequenceOverlaps().entrySet().stream()
                    .filter(x -> x.getKey() >= otherSeqIndex).forEach(x -> newSequenceOverlaps.put(x.getKey(), x.getValue()));
        }

        if(newSequence != null)
            primaryAssembly.updateSequenceInfo(newSequence, otherAssembly.sequenceOverlaps(), primaryOffsetAdjust);

        // merge read support, adjusting assembly indices and orientations as required
        if(primaryOffsetAdjust > 0)
        {
            for(JunctionAssembly assembly : primaryAssembly.assemblies())
            {
                for(SupportRead read : assembly.support())
                {
                    int newAssemblyIndex = read.fullAssemblyIndexStart() + primaryOffsetAdjust;
                    read.setFullAssemblyInfo(newAssemblyIndex, read.fullAssemblyOrientation());
                }
            }
        }

        if(otherOrientation.isReverse() || otherOffsetAdjust > 0)
        {
            // for non-reversed, reads need to be shifted by the amount that the primary starts earlier
            // for reversed reads, first reverse their coords based on just the non-primary, and then shift accordingly
            for(JunctionAssembly assembly : otherAssembly.assemblies())
            {
                for(SupportRead read : assembly.support())
                {
                    Orientation readOrientation;
                    int newAssemblyIndex;

                    if(otherOrientation.isReverse())
                    {
                        int origAssemblyIndex = read.fullAssemblyIndexStart();
                        newAssemblyIndex = otherLength - origAssemblyIndex - read.baseLength();

                        readOrientation = read.fullAssemblyOrientation().opposite();
                    }
                    else
                    {
                        newAssemblyIndex = read.fullAssemblyIndexStart();
                        readOrientation = read.fullAssemblyOrientation();
                    }

                    newAssemblyIndex += otherOffsetAdjust;

                    read.setFullAssemblyInfo(newAssemblyIndex, readOrientation);
                }
            }
        }
    }

    public static PhaseSet selectPrimaryPhaseSet(final PhaseSet first, final PhaseSet second)
    {
        // for now just select based on total support but will reconsider other factors just as length, mismatches
        int firstSupport = first.assemblies().stream().mapToInt(x -> x.supportCount()).sum();
        int secondSupport = second.assemblies().stream().mapToInt(x -> x.supportCount()).sum();

        return firstSupport >= secondSupport ? first : second;
    }
}