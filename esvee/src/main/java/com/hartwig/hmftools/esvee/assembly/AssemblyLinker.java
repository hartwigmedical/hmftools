package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PHASED_ASSEMBLY_MAX_TI;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_REF_SIDE_SOFT_CLIPS;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.INDEL;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionSequence;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

public final class AssemblyLinker
{
    public static AssemblyLink tryAssemblyFacing(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(!first.junction().Chromosome.equals(second.junction().Chromosome))
            return null;

        if(first.junction().Orient == second.junction().Orient)
            return null;

        JunctionAssembly lower = first.junction().Position < second.junction().Position ? first : second;
        JunctionAssembly upper = first == lower ? second : first;

        if(!lower.junction().isReverse())
            return null;

        int linkDistance = upper.junction().Position - lower.junction().Position;

        if(linkDistance < 0 || linkDistance > PHASED_ASSEMBLY_MAX_TI)
            return null;

        if(!first.refSideSoftClips().isEmpty() && !second.refSideSoftClips().isEmpty())
        {
            // cannot have ref aligned bases run past the other junction
            if(!refSideSoftClipMatchesJunction(lower, upper.junction().Position))
                return null;

            if(!refSideSoftClipMatchesJunction(upper, lower.junction().Position))
                return null;
        }
        else
        {
            // must share a junction read & mate in each
            boolean matched = false;

            for(SupportRead support : first.support())
            {
                if(!support.type().isSplitSupport())
                    continue;

                if(second.support()
                        .stream().filter(x -> x.type() == SupportType.JUNCTION)
                        .anyMatch(x -> x.matchesFragment(support)))
                {
                    matched = true;
                    break;
                }
            }

            if(!matched)
                return null;
        }

        return new AssemblyLink(lower, upper, LinkType.FACING, "", "");
    }

    private static boolean refSideSoftClipMatchesJunction(final JunctionAssembly assembly, int otherJunctionPosition)
    {
        if(assembly.refSideSoftClips().stream().noneMatch(x -> abs(x.Position - otherJunctionPosition) <= PROXIMATE_REF_SIDE_SOFT_CLIPS))
            return false;

        int refAlignedPosition = assembly.isForwardJunction() ? assembly.minAlignedPosition() : assembly.maxAlignedPosition();
        return abs(refAlignedPosition - otherJunctionPosition) <= PROXIMATE_REF_SIDE_SOFT_CLIPS;
    }

    public static AssemblyLink formLink(
            final JunctionAssembly first, final JunctionAssembly second, final JunctionSequence firstSeq, final JunctionSequence secondSeq,
            int firstIndexStart, int secondIndexStart, boolean isSecondary)
    {
        int firstJunctionOffset = firstSeq.junctionIndex() - firstIndexStart;
        int secondJunctionOffset = secondSeq.junctionIndex() - secondIndexStart;

        int firstJunctionIndexInSecond = -1;

        int junctionOffsetDiff = 0;

        if(!firstSeq.Reversed && !secondSeq.Reversed)
        {
            firstJunctionIndexInSecond = secondIndexStart + firstJunctionOffset;

            junctionOffsetDiff = secondJunctionOffset - firstJunctionOffset - 1;

        }
        else if(secondSeq.Reversed)
        {
            int firstJunctionIndexInSecondReversed = secondIndexStart + firstJunctionOffset;
            firstJunctionIndexInSecond = secondSeq.indexReverted(firstJunctionIndexInSecondReversed);

            junctionOffsetDiff = firstJunctionIndexInSecond - second.junctionIndex() - 1;
        }
        else
        {
            firstJunctionIndexInSecond = secondIndexStart + firstJunctionOffset;

            junctionOffsetDiff = second.junctionIndex() - firstJunctionIndexInSecond - 1;
        }

        String extraBases = "";

        if(junctionOffsetDiff != 0)
        {
            // inserted bases for any same-orientation break could be reverse-complemented at either breakend,
            // the convention is to do that for the higher breakend
            if(firstSeq.Reversed)
            {
                // FIXME: doesn't handle overlaps, only inserts
                extraBases = first.formSequence(
                        first.junctionIndex() - (junctionOffsetDiff), first.junctionIndex() - 1);
            }
            else
            {
                int insertStartIndex, insertEndIndex;

                // note the first sequence can use its normal junction coords since it is not reversed
                if(junctionOffsetDiff > 0)
                {
                    insertStartIndex = first.junctionIndex() + 1;
                    insertEndIndex = min(insertStartIndex + junctionOffsetDiff, firstSeq.FullSequence.length());
                }
                else
                {
                    insertStartIndex = first.junctionIndex() + junctionOffsetDiff + 1;
                    insertEndIndex = first.junctionIndex() + 1;
                }

                if(insertStartIndex < 0 || insertEndIndex <= insertStartIndex || insertEndIndex > firstSeq.BaseLength)
                {
                    if(isSecondary)
                        return null;

                    SV_LOGGER.debug("asm({} & {}) invalid insert/overlap indexRange({} -> {} junctOffsetDiff={} firstIndex={} secIndex={}) on firstSeq({}) secSeq({})",
                            first.junction().coords(), second.junction().coords(), insertStartIndex, insertEndIndex, junctionOffsetDiff,
                            firstIndexStart, secondIndexStart, firstSeq, secondSeq);

                    // attempt to correct to get some sequence  but needs investigation
                    if(insertStartIndex < 0)
                        insertStartIndex = 0;

                    if(insertEndIndex >= firstSeq.BaseLength)
                        insertEndIndex = firstSeq.BaseLength - 1;
                }

                if(insertEndIndex > insertStartIndex)
                {
                    extraBases = firstSeq.FullSequence.substring(insertStartIndex, insertEndIndex);
                }
            }
        }

        String insertedBases = junctionOffsetDiff > 0 ? extraBases : "";
        String overlapBases = junctionOffsetDiff < 0 ? extraBases : "";

        return new AssemblyLink(first, second, LinkType.SPLIT, insertedBases, overlapBases);
    }

    public static final int MATCH_SUBSEQUENCE_LENGTH = 10;

    public static AssemblyLink tryAssemblyIndel(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(!assembly1.indel() || !assembly2.indel())
            return null;

        // have already confirmed they share reads, so now just check the indel coords match
        if(!assembly1.indelCoords().matches(assembly2.indelCoords()))
            return null;

        String insertedBases = assembly1.indelCoords().insertedBases();

        return new AssemblyLink(assembly1, assembly2, INDEL, insertedBases, "");
    }

    public static AssemblyLink tryAssemblyOverlap(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        return tryAssemblyOverlap(assembly1, assembly2, false);
    }

    public static AssemblyLink tryAssemblyOverlap(final JunctionAssembly assembly1, final JunctionAssembly assembly2, boolean isSecondary)
    {
        JunctionAssembly first, second;
        JunctionSequence firstSeq, secondSeq;
        boolean firstReversed = false;
        boolean secondReversed = false;

        if(assembly1.junction().Orient != assembly2.junction().Orient)
        {
            first = assembly1.junction().isForward() ? assembly1 : assembly2;
            second = assembly1.junction().isReverse() ? assembly1 : assembly2;
        }
        else
        {
            first = assembly1;
            second = assembly2;

            if(assembly1.junction().isForward())
                secondReversed = true;
            else
                firstReversed = true;
        }

        firstSeq = new JunctionSequence(first, firstReversed);
        secondSeq = new JunctionSequence(second, secondReversed);

        // start with a simple comparison looking for the first sequence around its junction in the second
        String firstJunctionSequence = firstSeq.junctionSequence();
        int firstJunctionSeqLength = firstJunctionSequence.length();

        int firstSeqIndexInSecond = secondSeq.FullSequence.indexOf(firstJunctionSequence);

        if(firstSeqIndexInSecond >= 0)
        {
            // simple sequence match, so can form a link between these 2 assemblies
            return formLink(first, second, firstSeq, secondSeq, firstSeq.junctionSeqStartIndex(), firstSeqIndexInSecond, isSecondary);
        }

        if(isSecondary)
            return null;

        // take a smaller sections of the first's junction sequence and try to find their start index in the second sequence
        int juncSeqStartIndex = 0;
        List<int[]> alternativeIndexStarts = Lists.newArrayList();
        int subSeqIterations = (int)floor(firstJunctionSeqLength / MATCH_SUBSEQUENCE_LENGTH);
        for(int i = 0; i < subSeqIterations; ++i) // being the total junction sequence length (ie 100) divided by the subsequence length
        {
            juncSeqStartIndex = i * MATCH_SUBSEQUENCE_LENGTH;
            int juncSeqEndIndex = juncSeqStartIndex + MATCH_SUBSEQUENCE_LENGTH;

            if(juncSeqEndIndex >= firstJunctionSeqLength)
                break;

            String firstSubSequence = firstJunctionSequence.substring(juncSeqStartIndex, juncSeqStartIndex + MATCH_SUBSEQUENCE_LENGTH);

            int secondSubSeqIndex = secondSeq.FullSequence.indexOf(firstSubSequence);

            if(secondSubSeqIndex < 0)
                continue;

            alternativeIndexStarts.add(new int[] {juncSeqStartIndex, secondSubSeqIndex});

            secondSubSeqIndex = secondSeq.FullSequence.indexOf(firstSubSequence, secondSubSeqIndex + MATCH_SUBSEQUENCE_LENGTH);

            while(secondSubSeqIndex >= 0)
            {
                alternativeIndexStarts.add(new int[] {juncSeqStartIndex, secondSubSeqIndex});
                secondSubSeqIndex = secondSeq.FullSequence.indexOf(firstSubSequence, secondSubSeqIndex + MATCH_SUBSEQUENCE_LENGTH);
            }
        }

        // now perform a full junction sequence search in the second using the sequence matching logic
        Set<Integer> testedOffsets = Sets.newHashSet();

        int minOverlapLength = min(min(first.extensionLength(), second.extensionLength()), ASSEMBLY_LINK_OVERLAP_BASES);

        for(int[] indexStarts : alternativeIndexStarts)
        {
            // find a comparison range that falls within both sequence's index range around the junction
            int firstJuncSeqMatchIndex = indexStarts[0];
            int secondMatchIndex = indexStarts[1];

            int matchOffset = secondMatchIndex - firstJuncSeqMatchIndex;

            if(testedOffsets.contains(matchOffset))
                continue;

            testedOffsets.add(matchOffset);

            int secondIndexStart = secondMatchIndex - firstJuncSeqMatchIndex;
            int secondIndexEnd = secondIndexStart + firstJunctionSeqLength - 1;

            int firstJuncIndexStart = 0;
            int firstJuncIndexEnd = firstJunctionSeqLength - 1;

            if(secondIndexStart < 0)
            {
                firstJuncIndexStart += -(secondIndexStart);
                secondIndexStart = 0;
            }

            // discount this match if the implied end of the match in the second sequence is past its ref base end
            if(secondIndexEnd >= secondSeq.BaseLength)
                continue;

            int firstIndexStart = firstJuncIndexStart + firstSeq.junctionSeqStartIndex();
            int firstIndexEnd = min(firstJuncIndexEnd + firstSeq.junctionSeqStartIndex(), firstSeq.BaseLength - 1);

            if(secondIndexEnd - secondIndexStart + 1 < minOverlapLength || firstIndexEnd - firstIndexStart + 1 < minOverlapLength)
                continue;

            int mismatchCount = SequenceCompare.compareSequences(
                    firstSeq.bases(), firstSeq.baseQuals(), firstIndexStart, firstIndexEnd, firstSeq.repeatInfo(),
                    secondSeq.bases(), secondSeq.baseQuals(), secondIndexStart, secondIndexEnd, secondSeq.repeatInfo(),
                    PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchCount <= PRIMARY_ASSEMBLY_MERGE_MISMATCH)
            {
                // could hold out for a better match if there was more than one, but seem unlikely?
                return formLink(first, second, firstSeq, secondSeq, firstIndexStart, secondIndexStart, false);
            }
        }

        return null;
    }
}
