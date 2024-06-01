package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.setReadFlag;
import static com.hartwig.hmftools.esvee.TestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.esvee.TestUtils.DEFAULT_NM;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.alignment.AlignData;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;

import htsjdk.samtools.SAMFlag;

public class AssemblyTestUtils
{
    public static JunctionAssembly createAssembly(
            final String chromosome, final int junctionPosition, final Orientation junctionOrientation,
            final String assemblyBases, final int junctionIndex)
    {
        Junction junction = new Junction(chromosome, junctionPosition, junctionOrientation);

        int baseLength = assemblyBases.length();
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(baseLength);

        JunctionAssembly assembly = new JunctionAssembly(junction, assemblyBases.getBytes(), baseQuals, junctionIndex);
        assembly.buildRepeatInfo();
        return assembly;
    }

    public static AssemblyAlignment createAssemblyAlignment(
            final RefGenomeInterface refGenome, final String chrStart, int posStart, Orientation orientStart,
            final String chrEnd, int posEnd, Orientation orientEnd, final String insertedBases, final String overlapBases)
    {
        // first a basic exact match junction
        Junction junctionStart = new Junction(chrStart, posStart, orientStart);
        Junction junctionEnd = new Junction(chrEnd, posEnd, orientEnd);

        int refBaseLength = 100;
        int extBaseLength = 50;

        String firstRefBases, secondRefBases;

        if(orientStart.isForward())
            firstRefBases = refGenome.getBaseString(chrStart, posStart - refBaseLength + 1, posStart);
        else
            firstRefBases = refGenome.getBaseString(chrStart, posStart, posStart + refBaseLength - 1);

        if(orientEnd.isForward())
            secondRefBases = refGenome.getBaseString(chrEnd, posEnd - refBaseLength + 1, posEnd);
        else
            secondRefBases = refGenome.getBaseString(chrEnd, posEnd, posEnd + refBaseLength - 1);

        String firstExtBases = insertedBases;
        String secondExtBases = insertedBases;

        if(orientStart == orientEnd)
        {
            if(orientEnd.isForward())
                firstExtBases += secondRefBases.substring(secondRefBases.length() - extBaseLength);
            else
                firstExtBases += secondRefBases.substring(0, extBaseLength);

            firstExtBases = insertedBases + Nucleotides.reverseComplementBases(firstExtBases);

            if(orientStart.isForward())
                secondExtBases += firstRefBases.substring(firstRefBases.length() - extBaseLength);
            else
                secondExtBases += firstRefBases.substring(0, extBaseLength);

            secondExtBases = Nucleotides.reverseComplementBases(insertedBases) + Nucleotides.reverseComplementBases(secondExtBases);
        }
        else
        {
            if(orientEnd.isForward())
                firstExtBases += secondRefBases.substring(secondRefBases.length() - extBaseLength);
            else
                firstExtBases += secondRefBases.substring(0, extBaseLength);

            if(orientStart.isForward())
                secondExtBases += firstRefBases.substring(firstRefBases.length() - extBaseLength);
            else
                secondExtBases += firstRefBases.substring(0, extBaseLength);
        }

        String firstAssemblyBases, secondAssemblyBases;
        int firstJunctionIndex, secondJunctionIndex;

        if(orientStart.isForward())
        {
            firstAssemblyBases = firstRefBases + firstExtBases;
            firstJunctionIndex = firstRefBases.length() - 1;
        }
        else
        {
            firstAssemblyBases = firstExtBases + firstRefBases;
            firstJunctionIndex = firstExtBases.length();
        }

        if(orientEnd.isForward())
        {
            secondAssemblyBases = secondRefBases + secondExtBases;
            secondJunctionIndex = secondRefBases.length() - 1;
        }
        else
        {
            secondAssemblyBases = secondExtBases + secondRefBases;
            secondJunctionIndex = secondExtBases.length();
        }

        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(junctionStart, firstAssemblyBases.getBytes(), baseQuals, firstJunctionIndex);
        JunctionAssembly secondAssembly = new JunctionAssembly(junctionEnd, secondAssemblyBases.getBytes(), baseQuals, secondJunctionIndex);

        AssemblyLink link = new AssemblyLink(secondAssembly, firstAssembly, LinkType.SPLIT, insertedBases, overlapBases);

        return new AssemblyAlignment(0, link);
    }

    public static AlignData createAlignment(
            final String chromosome, int posStart, int posEnd, int segStart, int segEnd, final String cigar)
    {
        int score = segEnd - segStart + 1;
        return createAlignment(chromosome, posStart, posEnd, false, DEFAULT_MAP_QUAL, score, segStart, segEnd, cigar, "");
    }

    public static AlignData createAlignment(
            final String chromosome, int posStart, int posEnd, boolean reversed, int mapQual, int score, int segStart, int segEnd,
            final String cigar, final String xaTag)
    {
        String mdTag = "";

        int flags = 0;

        if(reversed)
            flags = setReadFlag(flags, SAMFlag.READ_REVERSE_STRAND);

        // note: as per BWA conventions 1 will be subtracted from the segment end
        return new AlignData(
                new ChrBaseRegion(chromosome, posStart, posEnd), segStart, segEnd, mapQual,score, flags, cigar, DEFAULT_NM, xaTag, mdTag);
    }
}
