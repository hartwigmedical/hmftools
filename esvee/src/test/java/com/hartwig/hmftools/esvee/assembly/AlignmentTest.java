package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.prep.TestUtils.setReadFlag;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.alignment.AlignData;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.alignment.BreakendBuilder;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;

import org.junit.Test;

import htsjdk.samtools.SAMFlag;

public class AlignmentTest
{
    private final MockRefGenome mRefGenome;
    private int mAssemblyId;

    private static final String CHR_1_REF_BASES = REF_BASES_RANDOM_100 + REF_BASES_200 + REF_BASES_RANDOM_100;
    private static final String CHR_2_REF_BASES = REF_BASES_RANDOM_100 + REF_BASES_200 + REF_BASES_RANDOM_100; // matching for now

    private static final int DEFAULT_MAP_QUAL = 60;
    private static final int DEFAULT_NM = 0;

    public AlignmentTest()
    {
        mRefGenome = new MockRefGenome();
        mAssemblyId = 0;

        mRefGenome.RefGenomeMap.put(CHR_1, CHR_1_REF_BASES);
        mRefGenome.RefGenomeMap.put(CHR_2, CHR_2_REF_BASES);
    }

    @Test
    public void testHomology()
    {


    }

    @Test
    public void testBasicSvTypes()
    {
        AssemblyAlignment assemblyAlignment = createAssemblyAlignment(
                CHR_1, 200, POS_ORIENT, CHR_2, 250, NEG_ORIENT, "");

        BreakendBuilder breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        AlignData alignment1 = createAlignment(CHR_1, 101, 200, 0, 149, "100M50S");
        AlignData alignment2 = createAlignment(CHR_2, 250, 349, 150, 299, "50S100M");

        List<AlignData> alignments = List.of(alignment1, alignment2);

        // breakendBuilder.formBreakends(alignments);

    }

    private static AlignData createAlignment(
            final String chromosome, int posStart, int posEnd, int segStart, int segEnd, final String cigar)
    {
        int score = segEnd - segStart + 1;
        return createAlignment(chromosome, posStart, posEnd, false, DEFAULT_MAP_QUAL, score, segStart, segEnd, cigar);
    }

    private static AlignData createAlignment(
            final String chromosome, int posStart, int posEnd, boolean reversed, int mapQual, int score, int segStart, int segEnd,
            final String cigar)
    {
        String xaTag = "";
        String mdTag = "";

        int flags = 0;

        if(reversed)
            flags = setReadFlag(flags, SAMFlag.READ_REVERSE_STRAND);

        return new AlignData(
                new ChrBaseRegion(chromosome, posStart, posEnd), segStart, segEnd, mapQual,score, flags, cigar, DEFAULT_NM, xaTag, mdTag);
    }

    private AssemblyAlignment createAssemblyAlignment(
            final String chrStart, int posStart, byte orientStart,
            final String chrEnd, int posEnd, byte orientEnd, final String insertedBases)
    {
        // first a basic exact match junction
        Junction junctionStart = new Junction(chrStart, posStart, orientStart);
        Junction junctionEnd = new Junction(chrEnd, posEnd, orientEnd);

        int refBaseLength = 100;
        int extBaseLength = 50;

        String firstRefBases, secondRefBases;

        if(orientStart == POS_ORIENT)
            firstRefBases = mRefGenome.getBaseString(chrStart, posStart - refBaseLength + 1, posStart);
        else
            firstRefBases = mRefGenome.getBaseString(chrStart, posStart, posStart + refBaseLength - 1);

        if(orientEnd == POS_ORIENT)
            secondRefBases = mRefGenome.getBaseString(chrEnd, posEnd - refBaseLength + 1, posEnd);
        else
            secondRefBases = mRefGenome.getBaseString(chrEnd, posEnd, posEnd + refBaseLength - 1);

        String firstExtBases = insertedBases;
        String secondExtBases = insertedBases;

        if(orientStart == orientEnd)
        {
            if(orientEnd == POS_ORIENT)
                firstExtBases = secondRefBases.substring(secondRefBases.length() - extBaseLength);
            else
                firstExtBases = secondRefBases.substring(0, extBaseLength);

            firstExtBases = insertedBases + Nucleotides.reverseComplementBases(firstExtBases);

            if(orientStart == POS_ORIENT)
                secondExtBases = firstRefBases.substring(firstRefBases.length() - extBaseLength);
            else
                secondExtBases = firstRefBases.substring(0, extBaseLength);

            secondExtBases = Nucleotides.reverseComplementBases(insertedBases) + Nucleotides.reverseComplementBases(secondExtBases);
        }
        else
        {
            if(orientEnd == POS_ORIENT)
                firstExtBases = secondRefBases.substring(secondRefBases.length() - extBaseLength);
            else
                firstExtBases = secondRefBases.substring(0, extBaseLength);

            if(orientStart == POS_ORIENT)
                secondExtBases = firstRefBases.substring(firstRefBases.length() - extBaseLength);
            else
                secondExtBases = firstRefBases.substring(0, extBaseLength);
        }

        String firstAssemblyBases, secondAssemblyBases;
        int firstJunctionIndex, secondJunctionIndex;

        if(orientStart == POS_ORIENT)
        {
            firstAssemblyBases = firstRefBases + firstExtBases;
            firstJunctionIndex = firstRefBases.length() - 1;
        }
        else
        {
            firstAssemblyBases = firstExtBases = firstRefBases;
            firstJunctionIndex = firstExtBases.length();
        }

        if(orientEnd == POS_ORIENT)
        {
            secondAssemblyBases = secondRefBases + secondExtBases;
            secondJunctionIndex = secondRefBases.length() - 1;
        }
        else
        {
            secondAssemblyBases = secondExtBases = secondRefBases;
            secondJunctionIndex = secondExtBases.length();
        }

        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(junctionStart, firstAssemblyBases.getBytes(), baseQuals, firstJunctionIndex);
        JunctionAssembly secondAssembly = new JunctionAssembly(junctionEnd, secondAssemblyBases.getBytes(), baseQuals, secondJunctionIndex);


        int firstJunctionIndexInSecond = 0;

        AssemblyLink link = new AssemblyLink(
                secondAssembly, firstAssembly, LinkType.SPLIT, firstJunctionIndexInSecond, insertedBases, "");

        return new AssemblyAlignment(mAssemblyId++, link);
    }
}
