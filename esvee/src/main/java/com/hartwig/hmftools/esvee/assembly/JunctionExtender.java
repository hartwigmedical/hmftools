package com.hartwig.hmftools.esvee.assembly;

import java.util.List;

import com.hartwig.hmftools.esvee.types.DiscordantGroup;
import com.hartwig.hmftools.esvee.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.types.JunctionExtension;
import com.hartwig.hmftools.esvee.read.Read;

public class JunctionExtender
{
    private final JunctionAssembly mAssembly;
    private final List<JunctionAssembly> mUnlinkedAssemblies;
    private final List<DiscordantGroup> mDiscordantGroups;
    private final JunctionExtension mExtension;

    public JunctionExtender(
            final JunctionAssembly assembly, final List<JunctionAssembly> unlinkedAssemblies, final List<DiscordantGroup> discordantGroups)
    {
        mAssembly = assembly;
        mDiscordantGroups = discordantGroups;
        mUnlinkedAssemblies = unlinkedAssemblies;
        mExtension = new JunctionExtension(mAssembly);
    }

    public void extendAssembly()
    {
        // find the read or assembly with the most overlap and base concordance, and use this to extend the sequence




        /*
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
         */


    }

    private boolean checkRead(final Read read)
    {


        return false;
    }

    private boolean checkDiscordantGroup(final DiscordantGroup discordantGroup)
    {


        return false;
    }


}
