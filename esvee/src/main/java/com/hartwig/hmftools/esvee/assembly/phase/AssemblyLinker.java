package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_LINK_DISC_ONLY_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PHASED_ASSEMBLY_MAX_TI;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PHASED_ASSEMBLY_MIN_TI;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_REF_SIDE_SOFT_CLIPS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.extractInsertSequence;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.tryLineSequenceLink;
import static com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder.isLocalAssemblyCandidate;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LOCAL_INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.JunctionSequence.PHASED_ASSEMBLY_MATCH_SEQ_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.SequenceCompare;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionSequence;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.common.IndelCoords;

public final class AssemblyLinker
{
    public static AssemblyLink tryAssemblyFacing(
            final JunctionAssembly first, final JunctionAssembly second, final List<AssemblyLink> splitLinks)
    {
        // cannot already be in a split link ie the DUP case
        if(splitLinks.stream().anyMatch(x -> x.matches(first, second)))
            return null;

        if(!first.junction().Chromosome.equals(second.junction().Chromosome))
            return null;

        if(first.junction().Orient == second.junction().Orient)
            return null;

        JunctionAssembly lower = first.junction().Position < second.junction().Position ? first : second;
        JunctionAssembly upper = first == lower ? second : first;

        if(!lower.junction().isReverse())
            return null;

        int linkDistance = upper.junction().Position - lower.junction().Position + 1;

        if(linkDistance < PHASED_ASSEMBLY_MIN_TI || linkDistance > PHASED_ASSEMBLY_MAX_TI)
            return null;

        boolean requireSharedRead = true;

        if(!first.refSideSoftClips().isEmpty() && !second.refSideSoftClips().isEmpty())
        {
            // cannot have ref aligned bases run past the other junction
            if(refSideSoftClipMatchesJunction(lower, upper.junction().Position)
            && refSideSoftClipMatchesJunction(upper, lower.junction().Position))
            {
                requireSharedRead = false;
            }
        }

        if(requireSharedRead)
        {
            if(first.indel() || second.indel())
            {
                // if they share a read and the read contains the indel coords, then consider this a facing link
                IndelCoords firstIndelCoords = first.indelCoords();
                IndelCoords secondIndelCoords = second.indelCoords();

                // must share a junction read and/or mate in each with one matching the indel coordinates
                boolean matched = false;

                for(SupportRead support : first.support())
                {
                    if(!support.type().isSplitSupport())
                        continue;

                    for(SupportRead secondSupport : second.support())
                    {
                        if(!secondSupport.type().isSplitSupport())
                            continue;

                        if(!secondSupport.matchesFragment(support, true))
                            continue;

                        if(firstIndelCoords != null)
                        {
                            if(support.indelCoords() != null && support.indelCoords().matches(firstIndelCoords))
                                matched = true;
                            else if(secondSupport.indelCoords() != null && secondSupport.indelCoords().matches(firstIndelCoords))
                                matched = true;
                        }
                        else
                        {
                            if(support.indelCoords() != null && support.indelCoords().matches(secondIndelCoords))
                                matched = true;
                            else if(secondSupport.indelCoords() != null && secondSupport.indelCoords().matches(secondIndelCoords))
                                matched = true;
                        }

                        if(matched)
                            break;
                    }
                }

                if(!matched)
                    return null;
            }
            else
            {
                // cannot have ref bases extending past each other's junctions
                if(lower.refBaseLength() > linkDistance || upper.refBaseLength() > linkDistance)
                    return null;

                boolean matched = false;

                // require a shared split, concordant read in either of the assemblies
                for(int i = 0; i <= 0; ++i)
                {
                    List<SupportRead> splitSupport = (i == 0) ? first.support() : second.support();
                    List<SupportRead> otherSupport = (i == 0) ? second.support() : first.support();

                    for(SupportRead support : splitSupport)
                    {
                        if(!support.type().isSplitSupport())
                            continue;

                        for(SupportRead other : otherSupport)
                        {
                            if(other.type() != JUNCTION && other.type() != DISCORDANT)
                                continue;

                            // must be a concordant pair unless the same read is both
                            if(!other.id().equals(support.id()))
                                continue;

                            if(other.flags() == support.flags()) // the same read
                            {
                                matched = true;
                            }
                            else if(other.orientation() != support.orientation())
                            {
                                // reads must face each other
                                if(other.orientation().isForward() && other.alignmentStart() < support.alignmentStart())
                                    matched = true;
                                else if(support.orientation().isForward() && support.alignmentStart() < other.alignmentStart())
                                    matched = true;
                            }

                            if(matched)
                                break;
                        }
                    }
                }

                if(!matched)
                    return null;
            }
        }

        // ensure the ref base positions of each assembly now match
        lower.trimRefBasePosition(upper.junction().Position);
        upper.trimRefBasePosition(lower.junction().Position);

        return new AssemblyLink(lower, upper, LinkType.FACING, "", "");
    }

    public static boolean isFacingAssemblyCandidate(
            final JunctionAssembly assembly, final Set<JunctionAssembly> facingAssemblies, final List<AssemblyLink> splitLinks)
    {
        if(facingAssemblies.contains(assembly))
            return false;

        if(assembly.outcome() == LOCAL_INDEL) // observed very few of these so excluded
            return false;

        if(assembly.discordantOnly())
        {
            // must have been linked
            if(splitLinks.stream().noneMatch(x -> x.hasAssembly(assembly)))
                return false;
        }

        return true;
    }

    private static boolean refSideSoftClipMatchesJunction(final JunctionAssembly assembly, int otherJunctionPosition)
    {
        if(assembly.refSideSoftClips().stream().noneMatch(x -> x.hasProximateMatch(otherJunctionPosition)))
            return false;

        return abs(assembly.refBasePosition() - otherJunctionPosition) <= PROXIMATE_REF_SIDE_SOFT_CLIPS;
    }

    public static boolean isAssemblyIndelLink(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        return assembly1.indel() && assembly2.indel()
            && assembly1.isForwardJunction() != assembly2.isForwardJunction()
            && assembly1.indelCoords().matches(assembly2.indelCoords());
    }

    public static AssemblyLink tryAssemblyIndel(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(!isAssemblyIndelLink(assembly1, assembly2))
            return null;

        String insertedBases = assembly1.indelCoords().insertedBases();

        return new AssemblyLink(assembly1, assembly2, INDEL, insertedBases, "");
    }

    @VisibleForTesting
    public static AssemblyLink tryAssemblyOverlap(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        boolean isLocalLink = isLocalAssemblyCandidate(assembly1, assembly2);
        return tryAssemblyOverlap(assembly1, assembly2, true, isLocalLink);
    }

    public static AssemblyLink tryAssemblyOverlap(
            final JunctionAssembly assembly1, final JunctionAssembly assembly2, boolean allowMismatches, boolean isLocalIndel)
    {
        JunctionAssembly first, second;
        JunctionSequence firstSeq, secondSeq;
        boolean firstReversed = false;
        boolean secondReversed = false;

        // select assemblies such that the sequence comparison is always comparing first in +ve orientation overlapping with
        // second in the -ve orientation, even one of them needs to be reverse complimented
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

        AssemblyLink lineLink = tryLineSequenceLink(first, second, firstReversed, secondReversed);

        if(lineLink != null)
            return lineLink;

        boolean requireRefBaseOverlap = isLocalIndel;

        if(requireRefBaseOverlap)
        {
            // assembly pairs in a local INDEL configuration haven't been checked for shared fragments, so it's necessary that they check
            // for an overlap in their ref bases, not just that they have a common insert sequence
            firstSeq = JunctionSequence.formStraddlingMatchSequence(first);
            secondSeq = JunctionSequence.formStraddlingMatchSequence(second);
        }
        else
        {
            firstSeq = JunctionSequence.formOuterExtensionMatchSequence(first, firstReversed);
            secondSeq = JunctionSequence.formOuterExtensionMatchSequence(second, secondReversed);
        }

        // start with a simple comparison looking for the first sequence around its junction in the second
        String firstMatchSequence = firstSeq.matchSequence();
        int firstMatchSeqLength = firstMatchSequence.length();

        int firstSeqIndexInSecond = secondSeq.FullSequence.indexOf(firstMatchSequence);

        if(firstSeqIndexInSecond >= 0)
        {
            // simple sequence match, so can form a link between these 2 assemblies
            return formLink(first, second, firstSeq, secondSeq, firstSeq.matchSeqStartIndex(), firstSeqIndexInSecond, 0);
        }

        if(!allowMismatches)
            return null;

        // extend to use full extension sequence from the first assembly for subsequence testing
        if(first.extensionLength() > 2 * PHASED_ASSEMBLY_MATCH_SEQ_LENGTH)
        {
            firstSeq = JunctionSequence.formFullExtensionMatchSequence(first, firstReversed);
            firstMatchSequence = firstSeq.matchSequence();
            firstMatchSeqLength = firstMatchSequence.length();
        }

        // take a smaller sections of the first's junction sequence and try to find their start index in the second sequence
        int matchSeqStartIndex = 0;
        List<int[]> alternativeIndexStarts = Lists.newArrayList();
        Set<String> testedSequences = Sets.newHashSet();
        int subSeqIterations = (int)floor(firstMatchSeqLength / MATCH_SUBSEQUENCE_LENGTH);
        for(int i = 0; i < subSeqIterations; ++i) // being the total junction sequence length (ie 100) divided by the subsequence length
        {
            matchSeqStartIndex = i * MATCH_SUBSEQUENCE_LENGTH;
            int matchSeqEndIndex = matchSeqStartIndex + MATCH_SUBSEQUENCE_LENGTH;

            if(matchSeqEndIndex >= firstMatchSeqLength + 1)
                break;

            String firstSubSequence = firstMatchSequence.substring(matchSeqStartIndex, matchSeqEndIndex);

            if(testedSequences.contains(firstSubSequence))
                continue;

            testedSequences.add(firstSubSequence);

            int secondSubSeqIndex = secondSeq.FullSequence.indexOf(firstSubSequence);

            if(secondSubSeqIndex < 0)
                continue;

            // check if the first match sequence (ie including all extension bases) fits/matches within the second
            // and if not whether there is sufficient overlap
            int impliedSecondSeqMatchSeqStart = secondSubSeqIndex - matchSeqStartIndex;
            int impliedSecondSeqMatchSeqEnd = impliedSecondSeqMatchSeqStart + firstMatchSeqLength - 1;

            if(impliedSecondSeqMatchSeqEnd >= secondSeq.BaseLength)
            {
                int secondMatchLength = secondSeq.BaseLength - impliedSecondSeqMatchSeqStart - 1;

                if(secondMatchLength < ASSEMBLY_LINK_OVERLAP_BASES)
                    continue; // still possible that a later subsequence will match earlier
            }

            alternativeIndexStarts.add(new int[] {matchSeqStartIndex, secondSubSeqIndex});

            secondSubSeqIndex = secondSeq.FullSequence.indexOf(firstSubSequence, secondSubSeqIndex + MATCH_SUBSEQUENCE_LENGTH);

            while(secondSubSeqIndex >= 0)
            {
                impliedSecondSeqMatchSeqStart = secondSubSeqIndex - matchSeqStartIndex;
                if(impliedSecondSeqMatchSeqStart + firstMatchSeqLength > secondSeq.BaseLength)
                    break;

                alternativeIndexStarts.add(new int[] {matchSeqStartIndex, secondSubSeqIndex});
                secondSubSeqIndex = secondSeq.FullSequence.indexOf(firstSubSequence, secondSubSeqIndex + MATCH_SUBSEQUENCE_LENGTH);
            }
        }

        if(alternativeIndexStarts.isEmpty())
            return null;

        // now perform a full junction sequence search in the second using the sequence matching logic
        int minOverlapLength = min(min(first.extensionLength(), second.extensionLength()), ASSEMBLY_LINK_OVERLAP_BASES);

        if(first.discordantOnly() && second.discordantOnly())
            minOverlapLength = min(minOverlapLength, ASSEMBLY_LINK_DISC_ONLY_OVERLAP_BASES);

        int[] topMatchIndices;

        if(requireRefBaseOverlap)
        {
            topMatchIndices = findLocalSequenceMatch(firstSeq, secondSeq, minOverlapLength, alternativeIndexStarts);
        }
        else
        {
            topMatchIndices = findBestSequenceMatch(firstSeq, secondSeq, minOverlapLength, alternativeIndexStarts);
        }

        if(topMatchIndices != null)
        {
            int firstIndexStart = topMatchIndices[0];
            int secondIndexStart = topMatchIndices[1];
            int firstPreJuncMismatchDiff = topMatchIndices[2];

            return formLink(first, second, firstSeq, secondSeq, firstIndexStart, secondIndexStart, firstPreJuncMismatchDiff);
        }

        return null;
    }

    public static int[] findLocalSequenceMatch(
            final JunctionSequence firstSeq, final JunctionSequence secondSeq, int minOverlapLength, final List<int[]> alternativeIndexStarts)
    {
        if(alternativeIndexStarts.isEmpty())
            return null;

        int topMatchLength = 0;
        double topMatchMismatchPenalty = 0;
        int[] topMatchIndices = null;

        Set<Integer> testedOffsets = Sets.newHashSet();

        // take each of the subsequence match locations, build out a longer sequence around it to include all extension bases for
        // each assembly, capped by the other's ref bases and then run the sequence-matching routine
        for(int[] indexStarts : alternativeIndexStarts)
        {
            int firstMatchSeqMatchIndex = indexStarts[0];

            int secondMatchIndex = indexStarts[1];

            int matchOffset = secondMatchIndex - firstMatchSeqMatchIndex;

            // skip testing a comparison anchored around the same offsets between the 2 sequences
            if(testedOffsets.contains(matchOffset))
                continue;

            testedOffsets.add(matchOffset);

            int firstMatchIndex = firstMatchSeqMatchIndex + firstSeq.matchSeqStartIndex();

            // extend in each direction to the end of the applicable extension bases
            int minLowerExtension = min(secondMatchIndex, firstMatchIndex);
            int maxLowerExtension = min(firstSeq.BaseLength - firstMatchIndex - 1, secondSeq.BaseLength - secondMatchIndex - 1);

            int firstIndexStart = firstMatchIndex - minLowerExtension;
            int firstIndexEnd = firstMatchIndex + maxLowerExtension;

            int secondIndexStart = secondMatchIndex - minLowerExtension;
            int secondIndexEnd = secondMatchIndex + maxLowerExtension;

            int overlapLength = firstIndexEnd - firstIndexStart + 1;

            if(overlapLength < minOverlapLength)
                continue;

            double mismatchPenalty = SequenceCompare.compareSequences(
                    firstSeq.bases(), firstSeq.baseQuals(), firstIndexStart, firstIndexEnd, firstSeq.repeatInfo(),
                    secondSeq.bases(), secondSeq.baseQuals(), secondIndexStart, secondIndexEnd, secondSeq.repeatInfo(),
                    PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchPenalty > PRIMARY_ASSEMBLY_MERGE_MISMATCH)
                continue;

            if(overlapLength > topMatchLength || (overlapLength == topMatchLength && mismatchPenalty < topMatchMismatchPenalty))
            {
                // check that the matching section covers each sequence's ref bases as well as the extension bases
                int minRefOverlapLength = MIN_VARIANT_LENGTH / 2;

                // note for INDELs, the first sequence is always positive orientation
                int firstRefBaseOverlap = firstSeq.junctionIndex() - firstIndexStart + 1;
                int firstExtBaseOverlap = firstIndexEnd - firstSeq.junctionIndex();

                int secondRefBaseOverlap = secondIndexEnd - secondSeq.junctionIndex() + 1;
                int secondExtBaseOverlap = secondSeq.junctionIndex() - secondIndexStart;

                if(firstRefBaseOverlap < minRefOverlapLength && secondRefBaseOverlap < minRefOverlapLength)
                    continue;

                if(firstExtBaseOverlap < minRefOverlapLength && secondExtBaseOverlap < minRefOverlapLength)
                    continue;
            }

            topMatchLength = overlapLength;
            topMatchIndices = new int[] {firstIndexStart, secondIndexStart, 0};
            topMatchMismatchPenalty = mismatchPenalty;
        }

        return topMatchIndices;
    }

    public static int[] findBestSequenceMatch(
            final JunctionSequence firstSeq, final JunctionSequence secondSeq, int minOverlapLength, final List<int[]> alternativeIndexStarts)
    {
        if(alternativeIndexStarts.isEmpty())
            return null;

        int topMatchLength = 0;
        double topMatchMismatchPenalty = 0;
        int[] topMatchIndices = null;

        Set<Integer> testedOffsets = Sets.newHashSet();

        int firstJunctionSeqLength = min(firstSeq.matchSequence().length(), PHASED_ASSEMBLY_MATCH_SEQ_LENGTH);

        // take each of the subsequence match locations, build out a longer sequence around it and check for a match
        // then return the longest of these
        for(int[] indexStarts : alternativeIndexStarts)
        {
            // find a comparison range that falls within both sequence's index range around the junction
            int firstMatchSeqMatchIndex = indexStarts[0];
            int secondMatchIndex = indexStarts[1];

            int matchOffset = secondMatchIndex - firstMatchSeqMatchIndex;

            // skip testing a comparison anchored around the same offsets between the 2 sequences
            if(testedOffsets.contains(matchOffset))
                continue;

            testedOffsets.add(matchOffset);

            int secondIndexStart = secondMatchIndex - firstMatchSeqMatchIndex;
            int secondIndexEnd = secondIndexStart + firstJunctionSeqLength - 1;

            int firstMatchIndexStart = 0;
            int firstMatchIndexEnd = firstJunctionSeqLength - 1;

            if(secondIndexStart < 0)
            {
                firstMatchIndexStart += -(secondIndexStart);
                firstMatchIndexEnd += -(secondIndexStart);
                secondIndexEnd += -(secondIndexStart);
                secondIndexStart = 0;
            }

            if(secondIndexEnd >= secondSeq.BaseLength)
            {
                // trim the match end if the extension bases exceed the length of the second sequence's ref bases, as can happen for short TIs
                int endReduction = secondIndexEnd - secondSeq.BaseLength + 1;
                secondIndexEnd -= endReduction;
                firstMatchIndexEnd -= endReduction;
            }

            int firstIndexStart = firstMatchIndexStart + firstSeq.matchSeqStartIndex();
            int firstIndexEnd = min(firstMatchIndexEnd + firstSeq.matchSeqStartIndex(), firstSeq.BaseLength - 1);

            int overlapLength = min(firstIndexEnd - firstIndexStart + 1, secondIndexEnd - secondIndexStart + 1);

            if(overlapLength < minOverlapLength)
                continue;

            double mismatchPenalty = SequenceCompare.compareSequences(
                    firstSeq.bases(), firstSeq.baseQuals(), firstIndexStart, firstIndexEnd, firstSeq.repeatInfo(),
                    secondSeq.bases(), secondSeq.baseQuals(), secondIndexStart, secondIndexEnd, secondSeq.repeatInfo(),
                    PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchPenalty > PRIMARY_ASSEMBLY_MERGE_MISMATCH)
                continue;

            // try to extend the matched sequence in both directions
            int[] matchExtensions = tryExtendMatchedSequence(
                    firstSeq.bases(), firstIndexStart, firstIndexEnd, secondSeq.bases(), secondIndexStart, secondIndexEnd);

            firstIndexStart -= matchExtensions[0];
            secondIndexStart -= matchExtensions[0];
            overlapLength += matchExtensions[0] + matchExtensions[1];

            if(overlapLength > topMatchLength || (overlapLength == topMatchLength && mismatchPenalty < topMatchMismatchPenalty))
            {
                topMatchLength = overlapLength;
                topMatchIndices = new int[] {firstIndexStart, secondIndexStart, 0};
                topMatchMismatchPenalty = mismatchPenalty;
            }
        }

        return topMatchIndices;
    }

    private static int[] tryExtendMatchedSequence(
            final byte[] firstBases, int firstIndexStart, int firstIndexEnd, final byte[] secondBases, int secondIndexStart, int secondIndexEnd)
    {
        int earlierMatchedBases = 0;
        int latterMatchedBases = 0;

        int firstIndex = firstIndexStart - 1;
        int secondIndex = secondIndexStart - 1;

        while(firstIndex >= 0 && secondIndex >= 0)
        {
            if(firstBases[firstIndex] != secondBases[secondIndex])
                break;

            ++earlierMatchedBases;
            --firstIndex;
            --secondIndex;
        }

        firstIndex = firstIndexEnd + 1;
        secondIndex = secondIndexEnd + 1;

        while(firstIndex < firstBases.length && secondIndex < secondBases.length)
        {
            if(firstBases[firstIndex] != secondBases[secondIndex])
                break;

            ++latterMatchedBases;
            ++firstIndex;
            ++secondIndex;
        }

        return new int[] {earlierMatchedBases, latterMatchedBases};
    }

    protected static AssemblyLink formLink(
            final JunctionAssembly first, final JunctionAssembly second, final JunctionSequence firstSeq, final JunctionSequence secondSeq,
            int firstIndexStart, int secondIndexStart, int firstMismatchDiff)
    {
        // translate a sequence match & overlap into an assembly link, capturing any overlapping or inserted bases
        int firstJunctionOffset = firstSeq.junctionIndex() - firstIndexStart - firstMismatchDiff;
        int secondJunctionOffset = secondSeq.junctionIndex() - secondIndexStart;

        int firstJunctionIndexInSecond = secondIndexStart + firstJunctionOffset;

        // determine whether the junctions align exactly (junctionOffsetDiff = 0), or has an overlap (junctionOffsetDiff < 0)
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

        String insertedBases = "";
        String overlapBases = "";

        if(junctionOffsetDiff > 0)
        {
            boolean secondReversed = first.isForwardJunction() == second.isForwardJunction();
            insertedBases = extractInsertSequence(first, false, second, secondReversed, junctionOffsetDiff);
        }
        else if(junctionOffsetDiff < 0)
        {
            int overlapLength = abs(junctionOffsetDiff);

            if(overlapLength >= first.refBaseLength() || overlapLength >= second.refBaseLength())
                return null;

            overlapBases = extractOverlapBases(first, overlapLength);

            if(overlapBases == null)
            {
                overlapBases = extractOverlapBases(second, junctionOffsetDiff);

                if(overlapBases == null)
                {
                    SV_LOGGER.debug("asm({} & {}) invalid insert/overlap junctOffsetDiff({} firstIndex={} secIndex={}) on firstSeq({}) secSeq({})",
                            first.junction().coords(), second.junction().coords(), junctionOffsetDiff,
                            firstIndexStart, secondIndexStart, firstSeq, secondSeq);

                    // failure to extract the required bases invalidates the link
                    return null;
                }
            }
        }

        return new AssemblyLink(first, second, LinkType.SPLIT, insertedBases, overlapBases);
    }

    private static String extractOverlapBases(final JunctionAssembly assembly, int overlapLength)
    {
        int extraBasesStartIndex, extraBasesEndIndex;

        // the first assembly is always positive orientation and so the extra bases can use the junction offset diff value directly
        // around its junction index - with the exception being for when both assemblies are -ve orientation, in which case the
        // first assembly has been reversed for sequence matching, and so its insert/overlap base capture must be switched

        // the extra bases captured are by convention always taken from the first assembly and not reverse-complimented
        if(assembly.isForwardJunction())
        {
            // an overlap, so go back from the junction
            extraBasesStartIndex = assembly.junctionIndex() - overlapLength + 1;
            extraBasesEndIndex = assembly.junctionIndex();
        }
        else
        {
            extraBasesStartIndex = assembly.junctionIndex();
            extraBasesEndIndex = assembly.junctionIndex() + overlapLength - 1;
        }

        if(extraBasesStartIndex >= 0 && extraBasesEndIndex >= extraBasesStartIndex && extraBasesEndIndex < assembly.baseLength())
        {
            return assembly.formSequence(extraBasesStartIndex, extraBasesEndIndex);
        }

        return null;
    }
}
