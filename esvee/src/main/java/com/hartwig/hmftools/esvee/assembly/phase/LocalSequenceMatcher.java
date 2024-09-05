package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createMinBaseQuals;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.findBestSequenceMatch;
import static com.hartwig.hmftools.esvee.assembly.types.JunctionSequence.PHASED_ASSEMBLY_MATCH_SEQ_LENGTH;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
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

        int chromosomeEnd = mRefGenome.getChromosomeLength(assembly.junction().Chromosome);

        if(localRegionEnd >= chromosomeEnd)
            return null;

        byte[] refGenomeBases = mRefGenome.getBases(assembly.junction().Chromosome, localRegionStart, localRegionEnd);

        byte[] refBaseQuals = createMinBaseQuals(refGenomeBases.length);

        JunctionSequence assemblySeq = JunctionSequence.formStraddlingMatchSequence(
                assembly, false, PHASED_ASSEMBLY_MATCH_SEQ_LENGTH / 2, -1);

        Orientation localRefOrientation = assembly.junction().Orient.opposite();
        JunctionSequence localRefSeq = new JunctionSequence(refGenomeBases, refBaseQuals, localRefOrientation, false);

        // start with a simple comparison looking for the first sequence around its junction in the second
        String assemblyExtBases = assembly.formJunctionSequence();
        int assemblyJunctionSeqLength = assemblyExtBases.length();

        // first a simple local match
        int assemblySeqIndexInRef = localRefSeq.FullSequence.indexOf(assemblyExtBases);

        if(assemblySeqIndexInRef >= 0)
        {
            // simple sequence match, so can form a link between these 2 assemblies
            return formLocalLink(assembly, localRegionStart, assemblySeqIndexInRef);
        }

        int matchSeqStartIndex = 0;
        List<int[]> alternativeIndexStarts = Lists.newArrayList();

        int subsequenceLength = MATCH_SUBSEQUENCE_LENGTH;
        int subSeqIterations = (int)floor(assemblyJunctionSeqLength / subsequenceLength);

        for(int i = 0; i < subSeqIterations; ++i)
        {
            matchSeqStartIndex = i * subsequenceLength;
            int matchSeqEndIndex = matchSeqStartIndex + subsequenceLength;

            if(matchSeqEndIndex >= assemblyJunctionSeqLength)
                break;

            String assemblySubSequence = assemblyExtBases.substring(matchSeqStartIndex, matchSeqStartIndex + subsequenceLength);

            int refSubSeqIndex = localRefSeq.FullSequence.indexOf(assemblySubSequence);

            if(refSubSeqIndex < 0)
                continue;

            alternativeIndexStarts.add(new int[] {matchSeqStartIndex, refSubSeqIndex});

            refSubSeqIndex = localRefSeq.FullSequence.indexOf(assemblySubSequence, refSubSeqIndex + subsequenceLength);

            while(refSubSeqIndex >= 0)
            {
                alternativeIndexStarts.add(new int[] {matchSeqStartIndex, refSubSeqIndex});
                refSubSeqIndex = localRefSeq.FullSequence.indexOf(assemblySubSequence, refSubSeqIndex + subsequenceLength);
            }
        }

        int minOverlapLength = min(assembly.extensionLength(), ASSEMBLY_LINK_OVERLAP_BASES);

        int[] topMatchIndices = findBestSequenceMatch(assemblySeq, localRefSeq, minOverlapLength, alternativeIndexStarts);

        if(topMatchIndices != null)
        {
            int secondIndexStart = topMatchIndices[1];

            // return formLocalLink(assembly, assemblySeq, firstIndexStart, secondIndexStart);
            return formLocalLink(assembly, localRegionStart, secondIndexStart);
        }

        return null;
    }

    private AssemblyLink formLocalLink(
            final JunctionAssembly assembly, final int localRegionStart, int localRefIndexStart)
    {
        // create a simple local assembly to contain this link info
        int localRefJunctionPos = localRegionStart + localRefIndexStart;

        int localRefJunctionIndex;
        Orientation localRefOrientation;

        int extensionBaseLength = assembly.extensionLength();

        // the local assembly will set the ref bases to the original assembly's extension bases, and take a fixed amount of its ref bases
        // for its extension bases
        int fixedAssemblyRefBaseLength = 50;
        String assemblyRefBases = assembly.formRefBaseSequence(fixedAssemblyRefBaseLength);
        String assemblyExtBases = assembly.formJunctionSequence();

        if(assembly.isForwardJunction())
        {
            // the local assembly has the opposite, so in this case negative orientation
            localRefOrientation = REVERSE;
            localRefJunctionIndex = assemblyRefBases.length();
        }
        else
        {
            localRefOrientation = FORWARD;

            // add the assembly's extension base length to get to the local of the junction in this local sequence
            localRefJunctionPos += extensionBaseLength - 1;
            localRefJunctionIndex = assemblyExtBases.length() - 1;
        }

        Junction localRefJunction = new Junction(assembly.junction().Chromosome, localRefJunctionPos, localRefOrientation);

        // trim the local ref sequence to 100 bases from the junction - no point in storing and writing 1000 bases

        String localAssemblySequence;

        if(localRefOrientation.isForward())
            localAssemblySequence = assemblyExtBases + assemblyRefBases;
        else
            localAssemblySequence = assemblyRefBases + assemblyExtBases;

        byte[] localAssemblyBases = localAssemblySequence.getBytes();
        byte[] localAssemblyQuals = createMinBaseQuals(localAssemblyBases.length);

        JunctionAssembly localRefAssembly = new JunctionAssembly(
                localRefJunction, localAssemblyBases, localAssemblyQuals, Lists.newArrayList(), Lists.newArrayList());

        localRefAssembly.setJunctionIndex(localRefJunctionIndex);

        if(assembly.isForwardJunction())
            return new AssemblyLink(assembly, localRefAssembly, LinkType.SPLIT, "", "");
        else
            return new AssemblyLink(localRefAssembly, assembly, LinkType.SPLIT, "", "");
    }
}
