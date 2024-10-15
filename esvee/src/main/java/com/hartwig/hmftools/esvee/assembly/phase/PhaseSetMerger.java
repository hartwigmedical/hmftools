package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createMinBaseQuals;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.types.JunctionSequence.PHASED_ASSEMBLY_MATCH_SEQ_LENGTH;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isLineInsertPair;

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

            if(!first.hasValidAssemblyAlignment() || first.isSecondaryLineLink())
                continue;

            for(int j = i + 1; j < phaseSets.size(); ++j)
            {
                PhaseSet second = phaseSets.get(j);

                if(!second.hasValidAssemblyAlignment() || second.isSecondaryLineLink())
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
            Orientation otherOrientation = (i == 0) ? FORWARD : REVERSE;

            String otherSequence = otherOrientation == FORWARD ?
                    otherFullSequence : Nucleotides.reverseComplementBases(otherFullSequence);

            int[] matchIndices = findMergeIndices(primaryAssembly, primarySequence, otherAssembly, otherSequence);

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

    private static int[] findMergeIndices(
            final AssemblyAlignment firstAssembly, final String firstSequence,
            final AssemblyAlignment secondAssembly, final String secondSequence)
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

            if(overlapLength < PHASED_ASSEMBLY_MATCH_SEQ_LENGTH)
                continue;

            int firstMatchIndexEnd = firstMatchIndexStart + overlapLength - 1;
            int secondMatchIndexEnd = secondMatchIndexStart + overlapLength - 1;

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

            if(mismatchCount > permittedMismatches)
                continue;

            if(!canMergePhaseSets(
                    firstAssembly, firstMatchIndexStart, firstMatchIndexEnd, secondAssembly, secondMatchIndexStart, secondMatchIndexEnd))
            {
                continue;
            }

            return new int[] { firstMatchIndexStart, secondMatchIndexStart };
        }

        return null;
    }

    private static boolean canMergePhaseSets(
            final AssemblyAlignment first, int firstOverlapStart, int firstOverlapEnd,
            final AssemblyAlignment second, final int secondOverlapStart, final int secondOverlapEnd)
    {
        if(matchSequenceOverlapsJunctions(first, firstOverlapStart, firstOverlapEnd))
            return true;

        if(matchSequenceOverlapsJunctions(second, secondOverlapStart, secondOverlapEnd))
            return true;

        for(JunctionAssembly firstAssembly : first.assemblies())
        {
            boolean firstIsLineSite = firstAssembly.hasLineSequence();

            for(JunctionAssembly secondAssembly : second.assemblies())
            {
                if((firstIsLineSite || secondAssembly.hasLineSequence()) && isLineInsertPair(firstAssembly, secondAssembly))
                    return true;
            }
        }

        return false;
    }

    private static boolean matchSequenceOverlapsJunctions(final AssemblyAlignment assemblyAlignment, int overlapStart, int overlapEnd)
    {
        return assemblyAlignment.linkIndices().stream().anyMatch(x -> positionWithin(x, overlapStart, overlapEnd));
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
            // the non-primary assembly starts earlier than the primary so extend at the start
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
}
