package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConstants.PHASED_ASSEMBLY_JUNCTION_OVERLAP;
import static com.hartwig.hmftools.esvee.SvConstants.PHASED_ASSEMBLY_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.LinkType;

public final class AssemblyLinker
{
    private static final int SUBSEQUENCE_LENGTH = 10;

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
                insertedBases = first.formSequence(first.junctionIndex() - impliedInsertedBaseLength, first.junctionIndex() - 1);
            }
            else
            {
                insertedBases =
                        firstSeq.FullSequence.substring(first.junctionIndex() + 1, first.junctionIndex() + 1 + impliedInsertedBaseLength);
            }
        }

        return new AssemblyLink(first, second, LinkType.SPLIT, firstJunctionIndexInSecond, insertedBases);
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

        // take a smaller sections of the first sequence and try to find their start index in the second sequence
        int subSequenceCount = firstSeq.comparisonLength() / SUBSEQUENCE_LENGTH;

        List<int[]> alternativeIndexStarts = Lists.newArrayList();

        for(int i = 0; i < subSequenceCount; ++i)
        {
            int firstSubSeqStartIndex = i * SUBSEQUENCE_LENGTH;

            // must be within of the bounds of what will cover the junction
            if(firstSubSeqStartIndex < firstSeq.compareSeqStartIndex())
                continue;

            if(firstSubSeqStartIndex > firstSeq.junctionIndex())
                break;

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

            int mismatchCount = SequenceCompare.compareSequences(
                    first.bases(), first.baseQuals(), firstIndexStart, firstIndexEnd, first.repeatInfo(),
                    second.bases(), second.baseQuals(), secondIndexStart, secondIndexEnd, second.repeatInfo(), PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchCount <= PRIMARY_ASSEMBLY_MERGE_MISMATCH)
            {
                // match found..
                return formLink(first, second, firstSeq, secondSeq, firstIndexStart, secondIndexStart);
            }
        }

        return null;
    }

    private static final int REF_BASE_LENGTH_CAP = 200;
    private static final int COMPARISON_RANGE = PHASED_ASSEMBLY_OVERLAP_BASES - PHASED_ASSEMBLY_JUNCTION_OVERLAP;

    private class AssemblySequence
    {
        public final boolean Reversed;
        public final String FullSequence;

        public final int ExtensionLength;
        public final int RefBaseLength; // may be capped
        public final int BaseLength;

        private final int mJunctionIndex;

        // indices for the min-overlap around the junction (eg +/- 30 bases)
        private final int mJunctionSeqIndexStart;
        private final int mJunctionSeqIndexEnd;

        // furthest back index into ref such that can compare 100 bases and find the min overlap - ie typically 70 bases
        private final int mCompareSeqIndexStart;
        private final int mCompareSeqIndexEnd;

        public AssemblySequence(final JunctionAssembly assembly, boolean reverseCompliment)
        {
            Reversed = reverseCompliment;
            RefBaseLength = min(assembly.refBaseLength(), REF_BASE_LENGTH_CAP);
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

    private class AssemblyIndexRange
    {
        private final JunctionAssembly mAssembly;
        public final int RefBaseLength;
        public final int ExtensionLength; // not including the junction index/position
        public final int IndexStart;
        public final int IndexEnd;

        public AssemblyIndexRange(final JunctionAssembly assembly)
        {
            mAssembly = assembly;
            int junctionIndex = assembly.junctionIndex();
            ExtensionLength = min(assembly.extensionLength(), PHASED_ASSEMBLY_OVERLAP_BASES - PHASED_ASSEMBLY_JUNCTION_OVERLAP);
            RefBaseLength = min(assembly.refBaseLength(), PHASED_ASSEMBLY_OVERLAP_BASES - PHASED_ASSEMBLY_JUNCTION_OVERLAP);

            if(assembly.isForwardJunction())
            {
                IndexStart = junctionIndex - RefBaseLength + 1;
                IndexEnd = junctionIndex + ExtensionLength;
            }
            else
            {
                IndexStart = junctionIndex - ExtensionLength;
                IndexEnd = junctionIndex + RefBaseLength - 1;
            }
        }

        public int[] junctionSequenceIndexRange()
        {
            int midpointLength = PHASED_ASSEMBLY_OVERLAP_BASES / 2;
            int junctionIndex = mAssembly.junctionIndex();
            int refLength = min(RefBaseLength, midpointLength);
            int extensionLength = min(ExtensionLength, midpointLength);

            if(mAssembly.isForwardJunction())
            {
                return new int[] {junctionIndex - refLength + 1, junctionIndex + extensionLength };
            }
            else
            {
                return new int[] {junctionIndex - extensionLength, junctionIndex + refLength - 1 };
            }
        }

        public int totalLength() { return ExtensionLength + RefBaseLength;}

        public String toString()
        {
            return format("range(%d-%d len=%d) ext(%d) ref(%d)",
                    IndexStart, IndexEnd, totalLength(), ExtensionLength, RefBaseLength);
        }
    }
}
