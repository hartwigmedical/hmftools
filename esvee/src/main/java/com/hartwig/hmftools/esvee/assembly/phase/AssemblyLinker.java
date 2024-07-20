package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PHASED_ASSEMBLY_MAX_TI;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_REF_SIDE_SOFT_CLIPS;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.INDEL;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.SequenceCompare;
import com.hartwig.hmftools.esvee.assembly.SequenceDiffInfo;
import com.hartwig.hmftools.esvee.assembly.SequenceDiffType;
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

        int linkDistance = upper.junction().Position - lower.junction().Position + 1;

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
            // cannot have ref bases extending past each other's junctions
            if(lower.refBaseLength() > linkDistance || upper.refBaseLength() > linkDistance)
                return null;

            // must share a junction read & mate in each
            boolean matched = false;

            for(SupportRead support : first.support())
            {
                if(!support.type().isSplitSupport())
                    continue;

                if(second.support()
                        .stream().filter(x -> x.type() == SupportType.JUNCTION)
                        .anyMatch(x -> x.matchesFragment(support, true)))
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

        return abs(assembly.refBasePosition() - otherJunctionPosition) <= PROXIMATE_REF_SIDE_SOFT_CLIPS;
    }

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
        return tryAssemblyOverlap(assembly1, assembly2, true);
    }

    public static AssemblyLink tryAssemblyOverlap(final JunctionAssembly assembly1, final JunctionAssembly assembly2, boolean allowMismatches)
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
            return formLink(first, second, firstSeq, secondSeq, firstSeq.junctionSeqStartIndex(), firstSeqIndexInSecond, 0);
        }

        if(!allowMismatches)
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
        int minOverlapLength = min(min(first.extensionLength(), second.extensionLength()), ASSEMBLY_LINK_OVERLAP_BASES);

        int[] topMatchIndices = findBestSequenceMatch(firstSeq, secondSeq, minOverlapLength, alternativeIndexStarts);

        if(topMatchIndices != null)
        {
            int firstIndexStart = topMatchIndices[0];
            int secondIndexStart = topMatchIndices[1];
            int firstPreJuncMismatchDiff = topMatchIndices[2];

            return formLink(first, second, firstSeq, secondSeq, firstIndexStart, secondIndexStart, firstPreJuncMismatchDiff);
        }

        return null;
    }

    public static int[] findBestSequenceMatch(
            final JunctionSequence firstSeq, final JunctionSequence secondSeq, int minOverlapLength, final List<int[]> alternativeIndexStarts)
    {
        int topMatchLength = 0;
        int topMatchMismatches = 0;
        int[] topMatchIndices = null;

        Set<Integer> testedOffsets = Sets.newHashSet();

        int firstJunctionSeqLength = firstSeq.junctionSequence().length();

        // take each of the subsequence match locations, build out a longer sequence around it and check for a match
        // then return the longest of these
        for(int[] indexStarts : alternativeIndexStarts)
        {
            // find a comparison range that falls within both sequence's index range around the junction
            int firstJuncSeqMatchIndex = indexStarts[0];
            int secondMatchIndex = indexStarts[1];

            int matchOffset = secondMatchIndex - firstJuncSeqMatchIndex;

            // skip testing a comparison anchored around the same offsets bewteen the 2 sequences
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

            int overlapLength = min(firstIndexEnd - firstIndexStart + 1, secondIndexEnd - secondIndexStart + 1);

            if(overlapLength < minOverlapLength)
                continue;

            int mismatchCount = SequenceCompare.compareSequences(
                    firstSeq.bases(), firstSeq.baseQuals(), firstIndexStart, firstIndexEnd, firstSeq.repeatInfo(),
                    secondSeq.bases(), secondSeq.baseQuals(), secondIndexStart, secondIndexEnd, secondSeq.repeatInfo(),
                    PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchCount > PRIMARY_ASSEMBLY_MERGE_MISMATCH)
                continue;

            if(overlapLength > topMatchLength)
            {
                topMatchLength = overlapLength;
                topMatchIndices = new int[] {firstIndexStart, secondIndexStart, 0};
                topMatchMismatches = mismatchCount;
            }
        }

        // if there were mismatches, check for a need to factor this into the distance between the first start index and its junction
        if(topMatchMismatches > 0)
        {
            int firstIndexStart = topMatchIndices[0];
            int firstIndexEnd = firstIndexStart + topMatchLength - 1;
            int secondIndexStart = topMatchIndices[1];
            int secondIndexEnd = secondIndexStart + topMatchLength - 1;

            List<SequenceDiffInfo> mismatchDiffs = SequenceCompare.getSequenceMismatchInfo(
                firstSeq.bases(), firstSeq.baseQuals(), firstIndexStart, firstIndexEnd, firstSeq.repeatInfo(),
                secondSeq.bases(), secondSeq.baseQuals(), secondIndexStart, secondIndexEnd, secondSeq.repeatInfo());

            int firstJuncIndex = firstSeq.junctionIndex();
            int firstJunctionRepeatDiff = 0;

            for(SequenceDiffInfo diffInfo : mismatchDiffs)
            {
                if(diffInfo.Type == SequenceDiffType.BASE)
                    continue;

                if(diffInfo.Index > firstJuncIndex)
                    break;

                firstJunctionRepeatDiff += diffInfo.BaseCount;
            }

            topMatchIndices[2] = firstJunctionRepeatDiff;
        }

        return topMatchIndices;
    }

    protected static AssemblyLink formLink(
            final JunctionAssembly first, final JunctionAssembly second, final JunctionSequence firstSeq, final JunctionSequence secondSeq,
            int firstIndexStart, int secondIndexStart, int firstMismatchDiff)
    {
        // translate a sequence match & overlap into an assembly link, capturing any overlapping or inserted bases
        int firstJunctionOffset = firstSeq.junctionIndex() - firstIndexStart - firstMismatchDiff;
        int secondJunctionOffset = secondSeq.junctionIndex() - secondIndexStart;

        int firstJunctionIndexInSecond = secondIndexStart + firstJunctionOffset;

        // determine whether the junctions align exactly (junctionOffsetDiff = 0), or have an overlap (junctionOffsetDiff < 0)
        // or there are inserted bases (junctionOffsetDiff > 0)
        int junctionOffsetDiff = 0;

        if(!firstSeq.Reversed && !secondSeq.Reversed)
        {
            junctionOffsetDiff = secondJunctionOffset - firstJunctionOffset - 1;

        }
        else if(secondSeq.Reversed)
        {
            firstJunctionIndexInSecond = secondSeq.indexReverted(firstJunctionIndexInSecond);
            junctionOffsetDiff = firstJunctionIndexInSecond - second.junctionIndex() - 1;
        }
        else
        {
            junctionOffsetDiff = second.junctionIndex() - firstJunctionIndexInSecond - 1;
        }

        String extraBases = "";

        if(junctionOffsetDiff != 0)
        {
            // the first assembly is always positive orientation and so the extra bases can use the junction offset diff value directly
            // around its junction index - with the exception being for when both assemblies are -ve orietation, in which case the
            // first assembly has been reversed for sequence matching, and so its insert/overlap base capture must be switched
            extraBases = extractExtraBases(first, firstSeq.Reversed, junctionOffsetDiff);

            if(extraBases == null)
            {
                extraBases = extractExtraBases(second, secondSeq.Reversed, junctionOffsetDiff);

                if(extraBases == null)
                {
                    SV_LOGGER.debug("asm({} & {}) invalid insert/overlap junctOffsetDiff({} firstIndex={} secIndex={}) on firstSeq({}) secSeq({})",
                            first.junction().coords(), second.junction().coords(), junctionOffsetDiff,
                            firstIndexStart, secondIndexStart, firstSeq, secondSeq);

                    extraBases = "";
                }
            }
        }

        String insertedBases = junctionOffsetDiff > 0 ? extraBases : "";
        String overlapBases = junctionOffsetDiff < 0 ? extraBases : "";

        return new AssemblyLink(first, second, LinkType.SPLIT, insertedBases, overlapBases);
    }

    private static String extractExtraBases(final JunctionAssembly assembly, boolean isReversed, int junctionOffsetDiff)
    {
        int extraBasesStartIndex, extraBasesEndIndex;

        // the first assembly is always positive orientation and so the extra bases can use the junction offset diff value directly
        // around its junction index - with the exception being for when both assemblies are -ve orietation, in which case the
        // first assembly has been reversed for sequence matching, and so its insert/overlap base capture must be switched

        // the extra bases captured are by convention always taken from the first assembly and not reverse-complimented
        if(assembly.isForwardJunction())
        {
            if(junctionOffsetDiff > 0)
            {
                // an insert, so take bases from after the junction
                extraBasesStartIndex = assembly.junctionIndex() + 1;
                extraBasesEndIndex = assembly.junctionIndex() + junctionOffsetDiff;
            }
            else
            {
                // an overlap, so go back from the junction
                extraBasesStartIndex = assembly.junctionIndex() + junctionOffsetDiff + 1;
                extraBasesEndIndex = assembly.junctionIndex();
            }
        }
        else
        {
            // now work around a -ve orientation junction
            if(junctionOffsetDiff > 0)
            {
                extraBasesStartIndex = assembly.junctionIndex() - junctionOffsetDiff;
                extraBasesEndIndex = assembly.junctionIndex() - 1;
            }
            else
            {
                extraBasesStartIndex = assembly.junctionIndex();
                extraBasesEndIndex = assembly.junctionIndex() - junctionOffsetDiff - 1;
            }
        }

        if(extraBasesStartIndex >= 0 && extraBasesEndIndex >= extraBasesStartIndex && extraBasesEndIndex < assembly.baseLength())
        {
            return assembly.formSequence(extraBasesStartIndex, extraBasesEndIndex);
        }

        return null;
    }
}
