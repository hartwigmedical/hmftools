package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.PHASED_ASSEMBLY_JUNCTION_OVERLAP;
import static com.hartwig.hmftools.esvee.SvConstants.PHASED_ASSEMBLY_MAX_TI;
import static com.hartwig.hmftools.esvee.SvConstants.PHASED_ASSEMBLY_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.SvConstants.PROXIMATE_REF_SIDE_SOFT_CLIPS;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findInsertedBases;
import static com.hartwig.hmftools.esvee.common.LinkType.INDEL;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.LinkType;
import com.hartwig.hmftools.esvee.common.RepeatInfo;
import com.hartwig.hmftools.esvee.read.Read;

public final class AssemblyLinker
{
    public static AssemblyLink tryAssemblyFacing(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(first.refSideSoftClips().isEmpty() || second.refSideSoftClips().isEmpty())
            return null;

        if(!first.junction().Chromosome.equals(second.junction().Chromosome))
            return null;

        if(first.junction().Orientation == second.junction().Orientation)
            return null;

        JunctionAssembly lower = first.junction().Position < second.junction().Position ? first : second;
        JunctionAssembly upper = first == lower ? second : first;

        if(!lower.junction().isReverse())
            return null;

        int linkDistance = upper.junction().Position - lower.junction().Position;

        if(linkDistance < 0 || linkDistance > PHASED_ASSEMBLY_MAX_TI)
            return null;

        // cannot have ref aligned bases run past the other junction
        if(!refSideSoftClipMatchesJunction(lower, upper.junction().Position))
            return null;

        if(!refSideSoftClipMatchesJunction(upper, lower.junction().Position))
            return null;

        return new AssemblyLink(lower, upper, LinkType.FACING, 0, "");
    }

    private static boolean refSideSoftClipMatchesJunction(final JunctionAssembly assembly, int otherJunctionPosition)
    {
        if(assembly.refSideSoftClips().stream().noneMatch(x -> abs(x.Position - otherJunctionPosition) <= PROXIMATE_REF_SIDE_SOFT_CLIPS))
            return false;

        int refAlignedPosition = assembly.isForwardJunction() ? assembly.minAlignedPosition() : assembly.maxAlignedPosition();
        return abs(refAlignedPosition - otherJunctionPosition) <= PROXIMATE_REF_SIDE_SOFT_CLIPS;
    }

    private static AssemblyLink formLink(
            final JunctionAssembly first, final JunctionAssembly second, final AssemblySequence firstSeq, final AssemblySequence secondSeq,
            int firstIndexStart, int secondIndexStart)
    {
        int firstJunctionOffset = firstSeq.junctionIndex() - firstIndexStart;
        int secondJunctionOffset = secondSeq.junctionIndex() - secondIndexStart;

        String insertedBases = "";
        int firstJunctionIndexInSecond = -1;
        int impliedInsertedBaseLength = 0;

        if(!firstSeq.Reversed && !secondSeq.Reversed)
        {
            firstJunctionIndexInSecond = secondIndexStart + firstJunctionOffset;

            impliedInsertedBaseLength = max(secondJunctionOffset - firstJunctionOffset - 1, 0);
        }
        else if(secondSeq.Reversed)
        {
            int firstJunctionIndexInSecondReversed = secondIndexStart + firstJunctionOffset;
            firstJunctionIndexInSecond = secondSeq.indexReverted(firstJunctionIndexInSecondReversed);

            impliedInsertedBaseLength = max(firstJunctionIndexInSecond - second.junctionIndex() - 1, 0);
        }
        else
        {
            firstJunctionIndexInSecond = secondIndexStart + firstJunctionOffset;

            impliedInsertedBaseLength = max(second.junctionIndex() - firstJunctionIndexInSecond - 1, 0);
        }

        if(impliedInsertedBaseLength > 0)
        {
            // inserted bases for any same-orientation break could be reverse-complemented at either breakend,
            // go with the convention of doing that for the higher breakend?
            if(firstSeq.Reversed)
            {
                insertedBases = first.formSequence(
                        first.junctionIndex() - impliedInsertedBaseLength, first.junctionIndex() - 1);
            }
            else
            {
                // TODO: any inserted bases should be retrievable if the junction index of one assembly has been found in the other
                // look into more examples from of the logs
                try
                {
                    int insertStartIndex = first.junctionIndex() + 1;
                    int insertEndIndex = min(insertStartIndex + impliedInsertedBaseLength, firstSeq.FullSequence.length());

                    if(insertEndIndex > insertStartIndex)
                    {
                        insertedBases = firstSeq.FullSequence.substring(insertStartIndex, insertEndIndex);
                    }
                }
                catch(Exception e)
                {
                    SV_LOGGER.debug("assembly({}) failed to form insert bases({}-{} len={}) firstJuncIndexInSec({}) from firstSeq({}) secondSeq({})",
                            first, first.junctionIndex() + 1, first.junctionIndex() + 1 + impliedInsertedBaseLength,
                            impliedInsertedBaseLength, firstJunctionIndexInSecond, firstSeq, secondSeq);
                }
            }
        }

        return new AssemblyLink(first, second, LinkType.SPLIT, firstJunctionIndexInSecond, insertedBases);
    }

    private static final int SUBSEQUENCE_LENGTH = 10;

    public AssemblyLink tryAssemblyIndel(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        if(!assembly1.indel() || !assembly2.indel())
            return null;

        // have already confirmed they share reads, so now just check the indel coords match
        if(!assembly1.initialRead().indelCoords().matches(assembly2.initialRead().indelCoords()))
            return null;

        String insertedBases = "";
        Read indelRead = assembly1.initialRead();

        if(indelRead.indelCoords().isInsert())
        {
            insertedBases = findInsertedBases(indelRead);
        }

        return new AssemblyLink(assembly1, assembly2, INDEL, 0, insertedBases);
    }

    public AssemblyLink tryAssemblyOverlap(final JunctionAssembly assembly1, final JunctionAssembly assembly2)
    {
        JunctionAssembly first, second;
        AssemblySequence firstSeq, secondSeq;
        boolean firstReversed = false;
        boolean secondReversed = false;

        if(assembly1.junction().Orientation != assembly2.junction().Orientation)
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

        firstSeq = new AssemblySequence(first, firstReversed);
        secondSeq = new AssemblySequence(second, secondReversed);

        // start with a simple comparison, taken from each of their junction index positions, ie assuming no inserted bases

        // first try a simple string search to find an overlap for the +/-30 bases around one assembly's junction in the other
        String firstJunctionSequence = firstSeq.junctionSequence();

        int firstSeqIndexInSecond = secondSeq.FullSequence.indexOf(firstJunctionSequence);

        if(firstSeqIndexInSecond >= 0)
        {
            // easy match, so can form a link at these 2 locations
            return formLink(first, second, firstSeq, secondSeq, firstSeq.junctionSeqStartIndex(), firstSeqIndexInSecond);
        }

        // AssemblyIndexRange secondIndexRange = new AssemblyIndexRange(second);

        // int subSequenceCount = firstSeq.comparisonLength() / SUBSEQUENCE_LENGTH;

        // take a smaller sections of the first sequence and try to find their start index in the second sequence
        List<int[]> alternativeIndexStarts = Lists.newArrayList();

        for(int firstSubSeqStartIndex = firstSeq.compareSeqStartIndex(); firstSubSeqStartIndex <= firstSeq.junctionIndex();
                firstSubSeqStartIndex += SUBSEQUENCE_LENGTH)
        {
            String firstSubSequence = firstSeq.FullSequence.substring(firstSubSeqStartIndex, firstSubSeqStartIndex + SUBSEQUENCE_LENGTH);

            int secondSubSeqIndex = secondSeq.FullSequence.indexOf(firstSubSequence);

            if(secondSubSeqIndex < 0)
                continue;

            // where must this match fall within the second's sequence? just need to ensure that the junction is covered in both
            if(!positionWithin(secondSubSeqIndex, secondSeq.compareSeqStartIndex(), secondSeq.compareSeqEndIndex()))
                continue;

            alternativeIndexStarts.add(new int[] {firstSubSeqStartIndex, secondSubSeqIndex});
        }

        // try each of these in turn for a full match - could take the one with the lowest if there are multiple?
        for(int[] indexStarts : alternativeIndexStarts)
        {
            // find a comparison range that falls within both sequence's index range around the junction
            int firstIndexStart = indexStarts[0];
            int secondIndexStart = indexStarts[1];

            if(secondIndexStart > second.junctionIndex())
            {
                // the comparison start points are shifted to cover the junctions in both
                int secondJuncOffset = second.junctionIndex() - secondIndexStart;
                int shiftDistance = min(abs(secondJuncOffset), firstIndexStart - firstSeq.compareSeqStartIndex());
                firstIndexStart -= shiftDistance;
                secondIndexStart -= shiftDistance;
            }

            int firstIndexEnd = firstIndexStart + PHASED_ASSEMBLY_OVERLAP_BASES - 1;
            int secondIndexEnd = secondIndexStart + PHASED_ASSEMBLY_OVERLAP_BASES - 1;

            if(firstIndexEnd >= first.baseLength() || secondIndexEnd >= second.baseLength())
                continue;

            int mismatchCount = SequenceCompare.compareSequences(
                    firstSeq.bases(), firstSeq.baseQuals(), firstIndexStart, firstIndexEnd, firstSeq.repeatInfo(),
                    secondSeq.bases(), secondSeq.baseQuals(), secondIndexStart, secondIndexEnd, secondSeq.repeatInfo(),
                    PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchCount <= PRIMARY_ASSEMBLY_MERGE_MISMATCH)
            {
                // could hold out for a better match if there was more than one, but seem unlikely?
                return formLink(first, second, firstSeq, secondSeq, firstIndexStart, secondIndexStart);
            }
        }

        return null;
    }

    private static final int REF_BASE_LENGTH_CAP = 100; // can't see why this would need to be any longer
    private static final int COMPARISON_RANGE = PHASED_ASSEMBLY_OVERLAP_BASES - PHASED_ASSEMBLY_JUNCTION_OVERLAP;

    private class AssemblySequence
    {
        public final boolean Reversed;
        public final String FullSequence;

        public final int ExtensionLength;
        public final int RefBaseLength; // may be capped
        public final int BaseLength;

        private final JunctionAssembly mAssembly;

        private final int mJunctionIndex;

        // indices for the min-overlap around the junction (eg +/- 30 bases)
        private final int mJunctionSeqIndexStart;
        private final int mJunctionSeqIndexEnd;

        // furthest back index into ref such that can compare 100 bases and find the min overlap - ie typically 70 bases
        private final int mCompareSeqIndexStart;
        private final int mCompareSeqIndexEnd;

        // built on demand since only used for the sequence comparison routine
        private List<RepeatInfo> mRepeatInfo;
        private byte[] mBases;
        private byte[] mBaseQuals;

        public AssemblySequence(final JunctionAssembly assembly, boolean reverseCompliment)
        {
            mAssembly = assembly;
            mBases = null;
            mBaseQuals = null;
            mRepeatInfo = null;

            Reversed = reverseCompliment;
            // RefBaseLength = min(assembly.refBaseLength(), REF_BASE_LENGTH_CAP); // not capped since throws out the junction index
            RefBaseLength = assembly.refBaseLength();
            ExtensionLength = assembly.extensionLength();
            BaseLength = RefBaseLength + ExtensionLength;

            mJunctionIndex = assembly.junctionIndex();

            if(!Reversed)
            {
                FullSequence = assembly.formJunctionSequence(RefBaseLength);
            }
            else
            {
                FullSequence = Nucleotides.reverseComplementBases(assembly.formJunctionSequence(RefBaseLength));
            }

            int compExtensionLength = min(ExtensionLength, COMPARISON_RANGE);
            int compRefBaseLength = min(RefBaseLength, COMPARISON_RANGE);

            int compIndexStart, compIndexEnd;

            if(assembly.isForwardJunction())
            {
                compIndexStart = mJunctionIndex - compRefBaseLength + 1;
                compIndexEnd = mJunctionIndex + compExtensionLength;
            }
            else
            {
                compIndexStart = mJunctionIndex - compExtensionLength;
                compIndexEnd = mJunctionIndex + compRefBaseLength - 1;
            }

            // also make a shorter sequence centred around the junction
            int juncSeqExtLength = min(ExtensionLength, PHASED_ASSEMBLY_JUNCTION_OVERLAP);
            int juncSeqRefLength = min(RefBaseLength, PHASED_ASSEMBLY_JUNCTION_OVERLAP);

            int juncIndexStart, juncIndexEnd;

            if(assembly.isForwardJunction())
            {
                juncIndexStart = mJunctionIndex - juncSeqRefLength + 1;
                juncIndexEnd = mJunctionIndex + juncSeqExtLength;
            }
            else
            {
                juncIndexStart = mJunctionIndex - juncSeqExtLength;
                juncIndexEnd = mJunctionIndex + juncSeqRefLength - 1;
            }

            if(!Reversed)
            {
                mJunctionSeqIndexStart = juncIndexStart;
                mJunctionSeqIndexEnd = juncIndexEnd;
                mCompareSeqIndexStart = compIndexStart;
                mCompareSeqIndexEnd = compIndexEnd;
            }
            else
            {
                // note the switches here
                mJunctionSeqIndexStart = indexReversed(juncIndexEnd);
                mJunctionSeqIndexEnd = indexReversed(juncIndexStart);
                mCompareSeqIndexStart = indexReversed(compIndexEnd);
                mCompareSeqIndexEnd = indexReversed(compIndexStart);
            }
        }

        public int junctionIndex() { return !Reversed ? mJunctionIndex : BaseLength - mJunctionIndex - 1; }

        public final String junctionSequence() { return FullSequence.substring(junctionSeqStartIndex(), junctionSeqEndIndex() + 1); }

        public int junctionSeqStartIndex() { return mJunctionSeqIndexStart; }
        public int junctionSeqEndIndex() { return mJunctionSeqIndexEnd; }
        public int compareSeqStartIndex() { return mCompareSeqIndexStart; }
        public int compareSeqEndIndex() { return mCompareSeqIndexEnd; }
        public int comparisonLength() { return mCompareSeqIndexEnd - mCompareSeqIndexStart + 1; }

        public byte[] bases()
        {
            if(!Reversed)
                return mAssembly.bases();

            if(mBases == null)
            {
                mBases = FullSequence.getBytes();
            }

            return mBases;
        }

        public byte[] baseQuals()
        {
            if(!Reversed)
                return mAssembly.baseQuals();

            if(mBaseQuals == null)
            {
                int baseLength = mAssembly.baseQuals().length;
                mBaseQuals = new byte[baseLength];

                for(int i = 0; i < baseLength; ++i)
                {
                    mBaseQuals[i] = mAssembly.baseQuals()[baseLength - i - 1];
                }
            }

            return mBaseQuals;
        }

        public List<RepeatInfo> repeatInfo()
        {
            if(mRepeatInfo == null)
            {
                List<RepeatInfo> repeats = RepeatInfo.findRepeats(FullSequence.getBytes());
                mRepeatInfo = repeats != null ? repeats : Collections.emptyList();
            }

            return mRepeatInfo;
        }

        private int indexReversed(int index)
        {
            int juncIndexDiff = index - mJunctionIndex;
            return junctionIndex() - juncIndexDiff;
        }

        public int indexReverted(int index)
        {
            int juncIndexDiff = index - junctionIndex();
            return mJunctionIndex - juncIndexDiff;
        }

        public String toString()
        {
            return format("len(%d ref=%d ext=%d juncIndex=%d) %s juncSeq(%d - %d) compSeq(%d - %d)",
                    BaseLength, RefBaseLength, ExtensionLength, mJunctionIndex, Reversed ? "rev" : "fwd",
                    junctionSeqStartIndex(), junctionSeqEndIndex(), compareSeqStartIndex(), compareSeqEndIndex());
        }
    }
}
