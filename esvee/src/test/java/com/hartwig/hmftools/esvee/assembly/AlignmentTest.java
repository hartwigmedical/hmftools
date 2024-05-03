package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.setReadFlag;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.createAssembly;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.alignment.HomologyData.determineHomology;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.alignment.AlignData;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.alignment.Breakend;
import com.hartwig.hmftools.esvee.alignment.BreakendBuilder;
import com.hartwig.hmftools.esvee.alignment.HomologyData;
import com.hartwig.hmftools.esvee.assembly.read.Read;
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
        AlignData alignmentStart = createAlignment(CHR_1, 100, 150, 0, 50, "51M");
        AlignData alignmentEnd = createAlignment(CHR_1, 100, 150, 51, 100, "50M");


        String assemblyOverlap = "";

        assertNull(determineHomology(assemblyOverlap, alignmentStart, alignmentEnd, mRefGenome));

        String basesStart = "TTCTTCTTCTC";
        String basesEnd = basesStart;
        assemblyOverlap = basesStart;

        // test 1: exact match
        HomologyData homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertNotNull(homology);
        assertEquals(basesStart, homology.Homology);
        assertEquals(-6, homology.ExactStart);
        assertEquals(5, homology.ExactEnd);
        assertEquals(-6, homology.InexactStart);
        assertEquals(5, homology.InexactEnd);

        // test 2: first base now no longer matches
        basesStart = "GTCTTCTTCTC";
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("", homology.Homology);
        assertEquals(0, homology.ExactStart);
        assertEquals(0, homology.ExactEnd);
        assertEquals(0, homology.InexactStart);
        assertEquals(11, homology.InexactEnd);

        // test 3: first base matches, range of lowest mismatches is 0-1
        basesStart = "TGCATCTTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = assemblyOverlap;
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("T", homology.Homology);
        assertEquals(-1, homology.ExactStart);
        assertEquals(0, homology.ExactEnd);
        assertEquals(-1, homology.InexactStart);
        assertEquals(10, homology.InexactEnd);

        // test 4: second base has mismatch
        basesStart = "TTGATCTTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = assemblyOverlap;
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("TT", homology.Homology);
        assertEquals(-1, homology.ExactStart);
        assertEquals(1, homology.ExactEnd);
        assertEquals(-1, homology.InexactStart);
        assertEquals(10, homology.InexactEnd);

        // test 5: min mismatches in range 4-6
        basesStart = "TTCATCGTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = "TTCTTCTTCTC";
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("TTC", homology.Homology);
        assertEquals(-1, homology.ExactStart);
        assertEquals(1, homology.ExactEnd);
        assertEquals(-5, homology.InexactStart);
        assertEquals(6, homology.InexactEnd);

        // test 6: ref sequences match but have the same difference from the assembly
        basesStart = "TTCATGTTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = basesStart;
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals(basesStart, homology.Homology);
        assertEquals(-6, homology.ExactStart);
        assertEquals(5, homology.ExactEnd);
        assertEquals(-6, homology.InexactStart);
        assertEquals(5, homology.InexactEnd);

        // test 7:
        basesStart = "TTCATGTTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = "TTCATATTCTC";
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("TTCAT", homology.Homology);
        assertEquals(-6, homology.ExactStart);
        assertEquals(5, homology.ExactEnd);
        assertEquals(-6, homology.InexactStart);
        assertEquals(5, homology.InexactEnd);

        // test 8: with an even number of bases
        basesStart = "AA";
        assemblyOverlap = basesStart;
        basesEnd = basesStart;
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals(basesStart, homology.Homology);
        assertEquals(-1, homology.ExactStart);
        assertEquals(1, homology.ExactEnd);
        assertEquals(-1, homology.InexactStart);
        assertEquals(1, homology.InexactEnd);
    }

    private static void addAssemblyRead(final JunctionAssembly assembly, int junctionOffset)
    {
        // a negative junction offset translates to soft-clip length
        int readStart;
        int softClip = junctionOffset < 0 ? abs(junctionOffset) : 0;
        String cigar;
        int readLength = 50;
        String readBases = REF_BASES_400.substring(0, readLength);

        if(assembly.isForwardJunction())
        {
            readStart = assembly.junction().Position - junctionOffset - readLength + 1;
            cigar = softClip > 0 ? format("%dM%dS", readLength - softClip, softClip) : format("%dM", readLength);
        }
        else
        {
            readStart = max(assembly.junction().Position + junctionOffset, assembly.junction().Position);
            cigar = softClip > 0 ? format("%dS%dM", softClip, readLength - softClip) : format("%dM", readLength);
        }

        Read read = createRead(READ_ID_GENERATOR.nextId(), readStart, readBases, cigar);
        assembly.addJunctionRead(read);
        assembly.setReadIndices();
    }

    @Test
    public void testAlignmentSequence()
    {
        String refBases1 = REF_BASES_400.substring(100, 200);
        String extBases1 = REF_BASES_400.substring(300, 350);
        String assemblyBases1 = refBases1 + extBases1;
        JunctionAssembly assembly1 = createAssembly(CHR_1, 200, FORWARD, assemblyBases1, refBases1.length() - 1);
        addAssemblyRead(assembly1, -30);
        addAssemblyRead(assembly1, 10);

        AssemblyAlignment assemblyAlignment = new AssemblyAlignment(0, assembly1);
        assertEquals(assemblyBases1, assemblyAlignment.fullSequence());
        assertEquals(150, assemblyAlignment.fullSequenceLength());
        assertEquals(80, assembly1.support().get(0).linkedAssemblyIndex());
        assertEquals(40, assembly1.support().get(1).linkedAssemblyIndex());

        refBases1 = REF_BASES_400.substring(100, 200);
        extBases1 = REF_BASES_400.substring(300, 350);
        assemblyBases1 = extBases1 + refBases1;
        assembly1 = createAssembly(CHR_1, 100, REVERSE, assemblyBases1, extBases1.length());

        assemblyAlignment = new AssemblyAlignment(0, assembly1);
        assertEquals(assemblyBases1, assemblyAlignment.fullSequence());
        assertEquals(150, assemblyAlignment.fullSequenceLength());

        // from a pair of assemblies
        refBases1 = REF_BASES_400.substring(101, 201);
        extBases1 = REF_BASES_400.substring(300, 350);
        assemblyBases1 = refBases1 + extBases1;
        assembly1 = createAssembly(CHR_1, 200, FORWARD, assemblyBases1, refBases1.length() - 1);

        String refBases2 = REF_BASES_400.substring(300, 400);
        String extBases2 = REF_BASES_400.substring(151, 201);
        String assemblyBases2 = extBases2 + refBases2;
        JunctionAssembly assembly2 = createAssembly(CHR_1, 300, REVERSE, assemblyBases2, extBases2.length());
        addAssemblyRead(assembly2, -30);
        assertEquals(20, assembly2.support().get(0).junctionAssemblyIndex());
        addAssemblyRead(assembly2, 10);
        assertEquals(60, assembly2.support().get(1).junctionAssemblyIndex());

        AssemblyLink assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, "", "");

        assemblyAlignment = new AssemblyAlignment(0, assemblyLink);
        String fullSequence = refBases1 + refBases2;
        assertEquals(fullSequence, assemblyAlignment.fullSequence());
        assertEquals(fullSequence.length(), assemblyAlignment.fullSequenceLength());
        assertEquals(70, assembly2.support().get(0).linkedAssemblyIndex());
        assertEquals(110, assembly2.support().get(1).linkedAssemblyIndex());

        // with overlap
        assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, "", refBases2.substring(0, 30));

        assemblyAlignment = new AssemblyAlignment(0, assemblyLink);
        fullSequence = refBases1 + refBases2.substring(30);
        assertEquals(fullSequence, assemblyAlignment.fullSequence());
        assertEquals(fullSequence.length(), assemblyAlignment.fullSequenceLength());

        // with inserted bases
        String insertedBases = "TTTTT";
        assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, insertedBases, "");

        assemblyAlignment = new AssemblyAlignment(0, assemblyLink);
        fullSequence = refBases1 + insertedBases + refBases2;
        assertEquals(fullSequence, assemblyAlignment.fullSequence());
        assertEquals(fullSequence.length(), assemblyAlignment.fullSequenceLength());

        // same again but with matching orientations
        assemblyBases1 = refBases1 + Nucleotides.reverseComplementBases(extBases1);
        assembly1 = createAssembly(CHR_1, 200, FORWARD, assemblyBases1, refBases1.length() - 1);

        assemblyBases2 = refBases2 + Nucleotides.reverseComplementBases(extBases2);
        assembly2 = createAssembly(CHR_1, 349, FORWARD, assemblyBases2, refBases2.length() - 1);

        assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, insertedBases, "");

        assemblyAlignment = new AssemblyAlignment(0, assemblyLink);
        fullSequence = refBases1 + insertedBases + Nucleotides.reverseComplementBases(refBases2);
        assertEquals(fullSequence, assemblyAlignment.fullSequence());
        assertEquals(fullSequence.length(), assemblyAlignment.fullSequenceLength());

        // both negative with an overlap
        String overlap = Nucleotides.reverseComplementBases(refBases2.substring(0, 10));
        extBases1 = Nucleotides.reverseComplementBases(refBases2.substring(10, 60));
        assemblyBases1 = extBases1 + refBases1;
        assembly1 = createAssembly(CHR_1, 101, REVERSE, assemblyBases1, extBases1.length());

        extBases2 = Nucleotides.reverseComplementBases(refBases1.substring(10, 60));
        assemblyBases2 = extBases2 + refBases2;
        assembly2 = createAssembly(CHR_1, 349, REVERSE, assemblyBases2, extBases2.length());

        assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, "", overlap);

        assemblyAlignment = new AssemblyAlignment(0, assemblyLink);
        fullSequence = Nucleotides.reverseComplementBases(refBases1) + refBases2.substring(10);
        assertEquals(fullSequence, assemblyAlignment.fullSequence());
        assertEquals(fullSequence.length(), assemblyAlignment.fullSequenceLength());
    }

    @Test
    public void testBasicSvTypes()
    {
        AssemblyAlignment assemblyAlignment = createAssemblyAlignment(
                CHR_1, 200, FORWARD, CHR_2, 250, REVERSE, "");

        BreakendBuilder breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        AlignData alignment1 = createAlignment(CHR_1, 101, 200, 0, 100, "100M50S");
        AlignData alignment2 = createAlignment(CHR_2, 250, 349, 100, 200, "50S100M");

        List<AlignData> alignments = Lists.newArrayList(alignment1, alignment2);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        Breakend first = assemblyAlignment.breakends().get(0);
        assertEquals(BND, first.svType());
        assertEquals(200, first.Position);
        assertEquals(FORWARD, first.Orient);
        assertNull(first.Homology);

        Breakend second = assemblyAlignment.breakends().get(1);
        assertEquals(250, second.Position);
        assertEquals(REVERSE, second.Orient);

        // test again for a DUP with inserted bases and SGLs on both ends
        assemblyAlignment = createAssemblyAlignment(
                CHR_1, 100, REVERSE, CHR_1, 350, FORWARD, "AAAA");

        breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        String altAlignment1 = "1,+1000000,100M,0;1,-1500000,50M50S,0";
        String altAlignment2 = "2,+2000000,100M,0";

        alignment1 = createAlignment(CHR_1, 251, 350, 0, 100, "100M50S");

        AlignData zeroAlignment1 = createAlignment(
                CHR_2, 1000, 1100, false, 0, 0, 100, 102, "100M", altAlignment1);

        AlignData zeroAlignment2 = createAlignment(
                CHR_2, 2000, 2100, false, 0, 0, 102, 104, "100M", altAlignment2);

        alignment2 = createAlignment(CHR_1, 100, 199, 104, 204, "50S100M");

        alignments = Lists.newArrayList(alignment1, alignment2, zeroAlignment1, zeroAlignment2);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        Collections.sort(assemblyAlignment.breakends());

        first = assemblyAlignment.breakends().get(0);
        assertEquals(DUP, first.svType());
        assertEquals(100, first.Position);
        assertEquals(REVERSE, first.Orient);
        assertEquals(5, first.alternativeAlignments().size());
        assertNull(first.Homology);

        second = assemblyAlignment.breakends().get(1);
        assertEquals(350, second.Position);
        assertEquals(FORWARD, second.Orient);
        assertEquals(5, second.alternativeAlignments().size());


        // test outer singles also with alt alignments
        assemblyAlignment = createAssemblyAlignment(
                CHR_1, 200, FORWARD, CHR_1, 250, REVERSE, "");

        breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        zeroAlignment1 = createAlignment(
                CHR_2, 3000, 3100, false, 0, 0, 0, 2, "50S2M", altAlignment1);

        alignment1 = createAlignment(CHR_1, 101, 200, 2, 101, "50S100M50S");

        alignment2 = createAlignment(CHR_1, 250, 349, 101, 199, "50S100M50S");

        zeroAlignment2 = createAlignment(
                CHR_2, 4000, 4100, false, 0, 0, 199, 200, "2M50S", altAlignment2);

        alignments = Lists.newArrayList(alignment1, alignment2, zeroAlignment1, zeroAlignment2);

        breakendBuilder.formBreakends(alignments);

        assertEquals(4, assemblyAlignment.breakends().size());

        Collections.sort(assemblyAlignment.breakends());

        Breakend breakend = assemblyAlignment.breakends().get(0);
        assertEquals(SGL, breakend.svType());
        assertEquals(101, breakend.Position);
        assertEquals(REVERSE, breakend.Orient);
        assertEquals(50, breakend.InsertedBases.length());
        assertEquals(3, breakend.alternativeAlignments().size());

        breakend = assemblyAlignment.breakends().get(1);
        assertEquals(DEL, breakend.svType());
        assertEquals(200, breakend.Position);
        assertEquals(FORWARD, breakend.Orient);
        assertNull(breakend.Homology);

        breakend = assemblyAlignment.breakends().get(2);
        assertEquals(DEL, breakend.svType());
        assertEquals(250, breakend.Position);
        assertEquals(REVERSE, breakend.Orient);

        breakend = assemblyAlignment.breakends().get(3);
        assertEquals(SGL, breakend.svType());
        assertEquals(349, breakend.Position);
        assertEquals(FORWARD, breakend.Orient);
        assertEquals(50, breakend.InsertedBases.length());
        assertEquals(2, breakend.alternativeAlignments().size());
    }

    @Test
    public void testIndelTypes()
    {
        // first a DEL
        AssemblyAlignment assemblyAlignment = createAssemblyAlignment(
                CHR_1, 200, FORWARD, CHR_1, 251, REVERSE, "");

        BreakendBuilder breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        AlignData alignment = createAlignment(CHR_1, 101, 350, 0, 200, "100M50D100M");

        List<AlignData> alignments = Lists.newArrayList(alignment);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        Breakend first = assemblyAlignment.breakends().get(0);
        assertEquals(DEL, first.svType());
        assertEquals(200, first.Position);
        assertEquals(FORWARD, first.Orient);
        assertNull(first.Homology);

        Breakend second = assemblyAlignment.breakends().get(1);
        assertEquals(251, second.Position);
        assertEquals(REVERSE, second.Orient);
    }

    private static AlignData createAlignment(
            final String chromosome, int posStart, int posEnd, int segStart, int segEnd, final String cigar)
    {
        int score = segEnd - segStart + 1;
        return createAlignment(chromosome, posStart, posEnd, false, DEFAULT_MAP_QUAL, score, segStart, segEnd, cigar, "");
    }

    private static AlignData createAlignment(
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

    public AssemblyAlignment createAssemblyAlignment(
            final String chrStart, int posStart, Orientation orientStart,
            final String chrEnd, int posEnd, Orientation orientEnd, final String insertedBases)
    {
        return createAssemblyAlignment(mRefGenome, chrStart, posStart, orientStart, chrEnd, posEnd, orientEnd, insertedBases);
    }

    public static AssemblyAlignment createAssemblyAlignment(
            final RefGenomeInterface refGenome, final String chrStart, int posStart, Orientation orientStart,
            final String chrEnd, int posEnd, Orientation orientEnd, final String insertedBases)
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

        AssemblyLink link = new AssemblyLink(secondAssembly, firstAssembly, LinkType.SPLIT, insertedBases, "");

        return new AssemblyAlignment(0, link);
    }
}
