package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConstants.PHASED_ASSEMBLY_JUNCTION_OVERLAP;
import static com.hartwig.hmftools.esvee.SvConstants.PHASED_ASSEMBLY_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.PrimaryPhaseGroup;

public class AssemblyOverlapper
{
    private final PrimaryPhaseGroup mPhaseGroup;

    public AssemblyOverlapper(final PrimaryPhaseGroup phaseGroup)
    {
        mPhaseGroup = phaseGroup;

    }

    public void buildPhasedAssembly()
    {
        List<JunctionAssembly> assemblies = mPhaseGroup.assemblies();

        if(assemblies.size() == 2)
        {
            tryAssemblyOverlap(assemblies.get(0), assemblies.get(1));
            return;
        }

        // where there are more than 2 assemblies, start with the ones with the most support and overlapping junction reads

        Collections.sort(assemblies, Comparator.comparingInt(x -> -x.supportCount()));

        List<JunctionAssembly> remainingAssemblies = assemblies.stream()
                .sorted(Comparator.comparingInt(x -> -x.supportCount()))
                .collect(Collectors.toList());

        // while(remainingAssemblies.isEmpty())


    }

    private class AssemblyIndexRange
    {
        private final JunctionAssembly mAssembly;
        public final int RefBaseLength;
        public final int ExtensionLength; // not including the junction index
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

    private static final int SUBSEQUENCE_LENGTH = 10;

    public boolean tryAssemblyOverlap(final JunctionAssembly first, final JunctionAssembly second)
    {
        AssemblyIndexRange firstIndexRange = new AssemblyIndexRange(first);

        // start with a simple comparison, taken from each of their junction index positions, ie assuming no inserted bases

        // for the case of matching orientations, one of the assemblies bases would need reversing - confirm this

        // first try a simple string search to find an overlap for the 100 bases around one assembly's junction in the other
        String firstFullSequence = first.formJunctionSequence(first.refBaseLength());

        int[] firstJunctionIndexRange = firstIndexRange.junctionSequenceIndexRange();
        String firstJunctionSequence = firstFullSequence.substring(firstJunctionIndexRange[0], firstJunctionIndexRange[1] + 1);

        String secondFullSequence = second.formJunctionSequence(second.refBaseLength());

        int firstIndexInSecond = secondFullSequence.indexOf(firstJunctionSequence);

        if(firstIndexInSecond >= 0)
        {
            // easy match, so can merge at these 2 locations...
            return true;
        }

        AssemblyIndexRange secondIndexRange = new AssemblyIndexRange(second);

        // take a smaller sections of the first sequence and try to find their start index in the second sequence
        int subSequenceCount = firstIndexRange.totalLength() / SUBSEQUENCE_LENGTH;

        List<int[]> alternativeIndexStarts = Lists.newArrayList();

        for(int i = 0; i < subSequenceCount; ++i)
        {
            int firstSubSeqStartIndex = i * SUBSEQUENCE_LENGTH;

            // must be within of the bounds of what will cover the junction
            if(firstSubSeqStartIndex < firstIndexRange.IndexStart)
                continue;

            //if(firstSubSeqStartIndex > first.junctionIndex() - PHASED_ASSEMBLY_JUNCTION_OVERLAP)
            if(firstSubSeqStartIndex > first.junctionIndex())
                break;

            String firstSubSequence = firstFullSequence.substring(firstSubSeqStartIndex, firstSubSeqStartIndex + SUBSEQUENCE_LENGTH);

            int secondSubSeqIndex = secondFullSequence.indexOf(firstSubSequence);

            if(secondSubSeqIndex < 0)
                continue;

            // where must this match fall within the second's sequence? just need to ensure that the junction is covered in both
            if(!positionWithin(secondSubSeqIndex, secondIndexRange.IndexStart, secondIndexRange.IndexEnd))
                continue;

            alternativeIndexStarts.add(new int[] {firstSubSeqStartIndex, secondSubSeqIndex});
        }

        // try each of these in turn for a full match
        for(int[] indexStarts : alternativeIndexStarts)
        {
            // find a comparison range that falls within both sequence's index range around the junction
            int firstIndexStart = indexStarts[0];
            int secondIndexStart = indexStarts[1];

            // int firstJuncOffset = first.junctionIndex() - firstIndexStart;

            if(secondIndexStart > second.junctionIndex())
            {
                // the comparison start points are shifted to cover the junctions in both
                int secondJuncOffset = second.junctionIndex() - secondIndexStart;
                int shiftDistance = min(abs(secondJuncOffset), firstIndexStart - firstIndexRange.IndexStart);
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
                return true;
            }
        }

        return false;
    }
}
