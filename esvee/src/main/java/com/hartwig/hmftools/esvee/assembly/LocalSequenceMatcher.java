package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.LOCAL_ASSEMBLY_REF_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PHASED_ASSEMBLY_JUNCTION_OVERLAP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyLinker.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyLinker.formLink;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createMinBaseQuals;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionSequence;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;

public class LocalSequenceMatcher
{
    private final RefGenomeInterface mRefGenome;
    private final int mLocalMatchDistance;

    public LocalSequenceMatcher(final RefGenomeInterface refGenome, final int localMatchDistance)
    {
        mRefGenome = refGenome;
        mLocalMatchDistance = localMatchDistance;
    }

    public AssemblyLink tryLocalAssemblyLink(final JunctionAssembly assembly)
    {
        int localRegionStart = assembly.junction().Position - mLocalMatchDistance;
        int localRegionEnd = assembly.junction().Position + mLocalMatchDistance;

        byte[] refGenomeBases = mRefGenome.getBases(assembly.junction().Chromosome, localRegionStart, localRegionEnd);

        byte[] refBaseQuals = createMinBaseQuals(refGenomeBases.length);

        JunctionSequence assemblySeq = new JunctionSequence(assembly, false, PHASED_ASSEMBLY_JUNCTION_OVERLAP, -1);

        byte localRefOrientation = assembly.isForwardJunction() ? NEG_ORIENT : POS_ORIENT;
        JunctionSequence localRefSeq = new JunctionSequence(refGenomeBases, refBaseQuals, localRefOrientation, false);

        // start with a simple comparison looking for the first sequence around its junction in the second
        String assemblyExtBases = assembly.formJunctionSequence();
        int assemblyJunctionSeqLength = assemblyExtBases.length();

        // first a simple local match
        int assemblySeqIndexInRef = localRefSeq.FullSequence.indexOf(assemblyExtBases);

        if(assemblySeqIndexInRef >= 0)
        {
            // simple sequence match, so can form a link between these 2 assemblies
            return formLocalLink(assembly, assemblySeq, localRegionStart, assemblySeqIndexInRef);
        }

        int juncSeqStartIndex = 0;
        List<int[]> alternativeIndexStarts = Lists.newArrayList();

        int subsequenceLength = MATCH_SUBSEQUENCE_LENGTH;
        int subSeqIterations = (int)floor(assemblyJunctionSeqLength / subsequenceLength);

        for(int i = 0; i < subSeqIterations; ++i)
        {
            juncSeqStartIndex = i * subsequenceLength;
            int juncSeqEndIndex = juncSeqStartIndex + subsequenceLength;

            if(juncSeqEndIndex >= assemblyJunctionSeqLength)
                break;

            String assemblySubSequence = assemblyExtBases.substring(juncSeqStartIndex, juncSeqStartIndex + subsequenceLength);

            int refSubSeqIndex = localRefSeq.FullSequence.indexOf(assemblySubSequence);

            if(refSubSeqIndex < 0)
                continue;

            alternativeIndexStarts.add(new int[] {juncSeqStartIndex, refSubSeqIndex});

            refSubSeqIndex = localRefSeq.FullSequence.indexOf(assemblySubSequence, refSubSeqIndex + subsequenceLength);

            while(refSubSeqIndex >= 0)
            {
                alternativeIndexStarts.add(new int[] {juncSeqStartIndex, refSubSeqIndex});
                refSubSeqIndex = localRefSeq.FullSequence.indexOf(assemblySubSequence, refSubSeqIndex + subsequenceLength);
            }
        }

        Set<Integer> testedOffsets = Sets.newHashSet();

        int minOverlapLength = min(assembly.extensionLength(), ASSEMBLY_LINK_OVERLAP_BASES);

        for(int[] indexStarts : alternativeIndexStarts)
        {
            // find a comparison range that falls within both sequence's index range around the junction
            int assemblyJuncSeqMatchIndex = indexStarts[0];
            int refMatchIndex = indexStarts[1];

            int matchOffset = refMatchIndex - assemblyJuncSeqMatchIndex;

            if(testedOffsets.contains(matchOffset))
                continue;

            testedOffsets.add(matchOffset);

            int refIndexStart = refMatchIndex - assemblyJuncSeqMatchIndex;
            int refIndexEnd = refIndexStart + assemblyJunctionSeqLength - 1;

            int assemblyJuncIndexStart = 0;
            int assemblyJuncIndexEnd = assemblyJunctionSeqLength - 1;

            if(refIndexStart < 0)
            {
                assemblyJuncIndexStart += -(refIndexStart);
                refIndexStart = 0;
            }

            // discount this match if the implied end of the match in the second sequence is past its ref base end
            if(refIndexEnd >= localRefSeq.BaseLength)
                continue;

            int assemblyIndexStart = assemblyJuncIndexStart + assemblySeq.junctionSeqStartIndex();
            int assemblyIndexEnd = min(assemblyJuncIndexEnd + assemblySeq.junctionSeqStartIndex(), assemblySeq.BaseLength - 1);

            if(refIndexEnd - refIndexStart + 1 < minOverlapLength || assemblyIndexEnd - assemblyIndexStart + 1 < minOverlapLength)
                continue;

            int mismatchCount = SequenceCompare.compareSequences(
                    assemblySeq.bases(), assemblySeq.baseQuals(), assemblyIndexStart, assemblyIndexEnd, assemblySeq.repeatInfo(),
                    localRefSeq.bases(), localRefSeq.baseQuals(), refIndexStart, refIndexEnd, localRefSeq.repeatInfo(),
                    PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchCount > PRIMARY_ASSEMBLY_MERGE_MISMATCH)
                continue;

            // now that the index in the remote ref sequence has a match and it is clear where this is in the assembly's extension sequence,
            // the implied junction position in the remote can be determined
            return formLocalLink(assembly, assemblySeq, localRegionStart, refIndexStart);
        }

        return null;
    }

    private AssemblyLink formLocalLink(
            final JunctionAssembly assembly, final JunctionSequence assemblySeq, final int localRegionStart, int localRefIndexStart)
    {
        // create a simple local assembly to contain this link info? depends perhaps on whether this will become an SV or not
        int localRefJunctionPos = localRegionStart + localRefIndexStart;

        int localRefSeqStart, localRefSeqEnd, localRefJunctionIndex;
        byte localRefOrientation;

        if(assembly.isForwardJunction())
        {
            // the local assembly has the opposite, so in this case negative orientation
            localRefOrientation = NEG_ORIENT;
            localRefSeqStart = localRefJunctionPos;
            localRefSeqEnd = localRefJunctionPos + LOCAL_ASSEMBLY_REF_LENGTH - 1;
            localRefJunctionIndex = 0;
        }
        else
        {
            localRefOrientation = POS_ORIENT;
            localRefSeqStart = localRefJunctionPos - LOCAL_ASSEMBLY_REF_LENGTH + 1;
            localRefSeqEnd = localRefJunctionPos;
            localRefJunctionIndex = LOCAL_ASSEMBLY_REF_LENGTH - 1;
        }

        Junction localRefJunction = new Junction(assembly.junction().Chromosome, localRefJunctionPos, localRefOrientation);

        // trim the local ref sequence to 100 bases from the junction - no point in storing and writing 1000 bases

        byte[] refGenomeBases = mRefGenome.getBases(assembly.junction().Chromosome, localRefSeqStart, localRefSeqEnd);
        byte[] refBaseQuals = createMinBaseQuals(refGenomeBases.length);

        JunctionAssembly localRefAssembly = new JunctionAssembly(
                localRefJunction, refGenomeBases, refBaseQuals, Lists.newArrayList(), Lists.newArrayList());

        localRefAssembly.setJunctionIndex(localRefJunctionIndex);

        return new AssemblyLink(assembly, localRefAssembly, LinkType.SPLIT, "", "");
    }
}
