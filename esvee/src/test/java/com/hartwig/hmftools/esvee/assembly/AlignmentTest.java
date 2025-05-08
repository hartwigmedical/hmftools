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
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.SSX2_GENE_ORIENT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.SSX2_MAX_MAP_QUAL;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.SSX2_REGIONS_V37;
import static com.hartwig.hmftools.esvee.TestUtils.DEFAULT_NM;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.alignment.AlignmentFilters.filterAlignments;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAlignment;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignData;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.alignment.Breakend;
import com.hartwig.hmftools.esvee.assembly.alignment.BreakendBuilder;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;

import org.junit.Test;

public class AlignmentTest
{
    private final MockRefGenome mRefGenome;

    private static final String CHR_1_REF_BASES = REF_BASES_RANDOM_100 + REF_BASES_400 + REF_BASES_RANDOM_100;
    private static final String CHR_2_REF_BASES = REF_BASES_RANDOM_100 + REF_BASES_400 + REF_BASES_RANDOM_100; // matches chr1
    private static final String CHR_3_REF_BASES = REF_BASES_RANDOM_100 + REF_BASES_200 + REF_BASES_RANDOM_100; // has much more repetition

    public AlignmentTest()
    {
        mRefGenome = new MockRefGenome();
        mRefGenome.RefGenomeMap.put(CHR_1, CHR_1_REF_BASES);
        mRefGenome.RefGenomeMap.put(CHR_2, CHR_2_REF_BASES);
        mRefGenome.RefGenomeMap.put(CHR_3, CHR_3_REF_BASES);
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

        AssemblyAlignment assemblyAlignment = new AssemblyAlignment(assembly1);
        assertEquals(assemblyBases1, assemblyAlignment.fullSequence());
        assertEquals(150, assemblyAlignment.fullSequenceLength());

        assertEquals(80, assembly1.support().get(0).fullAssemblyIndexStart());
        assertEquals(40, assembly1.support().get(1).fullAssemblyIndexStart());

        refBases1 = REF_BASES_400.substring(100, 200);
        extBases1 = REF_BASES_400.substring(300, 350);
        assemblyBases1 = extBases1 + refBases1;
        assembly1 = createAssembly(CHR_1, 100, REVERSE, assemblyBases1, extBases1.length());

        assemblyAlignment = new AssemblyAlignment(assembly1);
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
        assertEquals(30, assembly2.support().get(0).junctionReadStartDistance());
        addAssemblyRead(assembly2, 10);
        assertEquals(-10, assembly2.support().get(1).junctionReadStartDistance());

        AssemblyLink assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, "", "");

        assemblyAlignment = AssemblyTestUtils.createAssemblyAlignment(assemblyLink);
        String fullSequence = refBases1 + refBases2;
        assertEquals(fullSequence, assemblyAlignment.fullSequence());
        assertEquals(fullSequence.length(), assemblyAlignment.fullSequenceLength());
        assertEquals(70, assembly2.support().get(0).fullAssemblyIndexStart());
        assertEquals(110, assembly2.support().get(1).fullAssemblyIndexStart());

        // with overlap
        assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, "", refBases2.substring(0, 30));

        assemblyAlignment = AssemblyTestUtils.createAssemblyAlignment(assemblyLink);
        fullSequence = refBases1 + refBases2.substring(30);
        assertEquals(fullSequence, assemblyAlignment.fullSequence());
        assertEquals(fullSequence.length(), assemblyAlignment.fullSequenceLength());

        // with inserted bases
        String insertedBases = "TTTTT";

        assemblyBases1 = refBases1 + insertedBases + extBases1;
        assembly1 = createAssembly(CHR_1, 200, FORWARD, assemblyBases1, refBases1.length() - 1);
        assemblyBases2 = extBases2 + insertedBases + refBases2;
        assembly2 = createAssembly(CHR_1, 300, REVERSE, assemblyBases2, extBases2.length() + insertedBases.length());

        assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, insertedBases, "");

        assemblyAlignment = AssemblyTestUtils.createAssemblyAlignment(assemblyLink);
        fullSequence = refBases1 + insertedBases + refBases2;
        assertEquals(fullSequence, assemblyAlignment.fullSequence());
        assertEquals(fullSequence.length(), assemblyAlignment.fullSequenceLength());

        // same again but with matching orientations
        assemblyBases1 = refBases1 + insertedBases + Nucleotides.reverseComplementBases(extBases1);
        assembly1 = createAssembly(CHR_1, 200, FORWARD, assemblyBases1, refBases1.length() - 1);

        assemblyBases2 = refBases2 + Nucleotides.reverseComplementBases(insertedBases) + Nucleotides.reverseComplementBases(extBases2);
        assembly2 = createAssembly(CHR_1, 349, FORWARD, assemblyBases2, refBases2.length() - 1);

        assemblyLink = new AssemblyLink(assembly1, assembly2, LinkType.SPLIT, insertedBases, "");

        assemblyAlignment = AssemblyTestUtils.createAssemblyAlignment(assemblyLink);
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

        assemblyAlignment = AssemblyTestUtils.createAssemblyAlignment(assemblyLink);
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
        assertFalse(first.Homology.exists());

        Breakend second = assemblyAlignment.breakends().get(1);
        assertEquals(250, second.Position);
        assertEquals(REVERSE, second.Orient);

        // test again for a DUP with inserted bases
        assemblyAlignment = createAssemblyAlignment(
                CHR_1, 100, REVERSE, CHR_1, 350, FORWARD, "AAAA");

        breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        String altAlignment1 = "1,+1000000,100M,0;1,-1500000,50M50S,0";
        String altAlignment2 = "2,+2000000,100M,0";

        alignment1 = createAlignment(CHR_1, 251, 350, 0, 100, "100M50S");

        AlignData zeroAlignment1 = createAlignment(
                CHR_2, 200000, 200100, false, 0, 0, 100, 102, "100M",
                altAlignment1, "");

        AlignData zeroAlignment2 = createAlignment(
                CHR_2, 300000, 300100, false, 0, 0, 102, 104, "100M",
                altAlignment2, "");

        alignment2 = createAlignment(CHR_1, 100, 199, 104, 204, "50S100M");

        alignments = Lists.newArrayList(alignment1, alignment2);
        // alignments = Lists.newArrayList(alignment1, zeroAlignment1, zeroAlignment2, alignment2);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        Collections.sort(assemblyAlignment.breakends());

        first = assemblyAlignment.breakends().get(0);
        assertEquals(DUP, first.svType());
        assertEquals(100, first.Position);
        assertEquals(REVERSE, first.Orient);
        assertFalse(first.Homology.exists());

        second = assemblyAlignment.breakends().get(1);
        assertEquals(350, second.Position);
        assertEquals(FORWARD, second.Orient);
        // assertEquals(5, second.alternativeAlignments().size());


        // test outer singles also with alt alignments
        assemblyAlignment = createAssemblyAlignment(
                CHR_1, 200, FORWARD, CHR_1, 250, REVERSE, "");

        breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        zeroAlignment1 = createAlignment(
                CHR_2, 20000, 20100, false, 0, 0, 0, 2, "50S2M",
                altAlignment1, "");

        alignment1 = createAlignment(CHR_1, 101, 200, 2, 101, "50S100M50S");

        alignment2 = createAlignment(CHR_1, 250, 349, 101, 199, "50S100M50S");

        zeroAlignment2 = createAlignment(
                CHR_2, 40000, 40100, false, 0, 0, 199, 200, "2M50S",
                altAlignment2, "");

        alignments = Lists.newArrayList(zeroAlignment1, alignment1, alignment2, zeroAlignment2);
        alignments = Lists.newArrayList(alignment1, alignment2);

        breakendBuilder.formBreakends(alignments);

        assertEquals(4, assemblyAlignment.breakends().size());

        Collections.sort(assemblyAlignment.breakends());

        Breakend breakend = assemblyAlignment.breakends().get(0);
        assertEquals(SGL, breakend.svType());
        assertEquals(101, breakend.Position);
        assertEquals(REVERSE, breakend.Orient);
        assertEquals(50, breakend.InsertedBases.length());
        // assertEquals(3, breakend.alternativeAlignments().size());

        breakend = assemblyAlignment.breakends().get(1);
        assertEquals(DEL, breakend.svType());
        assertEquals(200, breakend.Position);
        assertEquals(FORWARD, breakend.Orient);
        assertFalse(first.Homology.exists());

        breakend = assemblyAlignment.breakends().get(2);
        assertEquals(DEL, breakend.svType());
        assertEquals(250, breakend.Position);
        assertEquals(REVERSE, breakend.Orient);

        breakend = assemblyAlignment.breakends().get(3);
        assertEquals(SGL, breakend.svType());
        assertEquals(349, breakend.Position);
        assertEquals(FORWARD, breakend.Orient);
        assertEquals(50, breakend.InsertedBases.length());
        // assertEquals(2, breakend.alternativeAlignments().size());
    }

    @Test
    public void testIndelTypes()
    {
        // first a DEL
        AssemblyAlignment assemblyAlignment = createAssemblyAlignment(
                CHR_3, 200, FORWARD, CHR_3, 251, REVERSE, "");

        BreakendBuilder breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        AlignData alignment = createAlignment(CHR_3, 101, 350, 0, 200, "100M50D100M");

        List<AlignData> alignments = Lists.newArrayList(alignment);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        Breakend first = assemblyAlignment.breakends().get(0);
        assertEquals(DEL, first.svType());
        assertEquals(200, first.Position);
        assertEquals(FORWARD, first.Orient);
        assertFalse(first.Homology.exists());

        Breakend second = assemblyAlignment.breakends().get(1);
        assertEquals(251, second.Position);
        assertEquals(REVERSE, second.Orient);


        // a DEL with homology
        assemblyAlignment = createAssemblyAlignment(
                CHR_3, 132, FORWARD, CHR_1, 190, REVERSE, "");

        breakendBuilder = new BreakendBuilder(mRefGenome, assemblyAlignment);

        // homology at pos 13 and 40
        alignment = createAlignment(CHR_3, 80, 250, 0, 120, "53M57D60M");

        alignments = Lists.newArrayList(alignment);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        first = assemblyAlignment.breakends().get(0);
        assertEquals(DEL, first.svType());
        assertEquals(134, first.Position);
        assertEquals(FORWARD, first.Orient);
        assertNotNull(first.Homology);
        assertEquals("CCCC", first.Homology.Homology);
        assertEquals(-2, first.Homology.ExactStart);
        assertEquals(2, first.Homology.ExactEnd);

        second = assemblyAlignment.breakends().get(1);
        assertEquals(192, second.Position);
        assertEquals(REVERSE, second.Orient);
    }

    @Test
    public void testAlignmentQualCalcs()
    {
        AssemblyAlignment assemblyAlignment = createAssemblyAlignment(
                CHR_1, 200, FORWARD, CHR_2, 200, REVERSE, "");

        AlignData alignment1 = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 200), 0, 100,
                60, 100, 0, "100M", DEFAULT_NM, "", "");

        AlignData alignment2 = new AlignData(
                new ChrBaseRegion(CHR_2, 1, 100), 100, 200,
                60, 100, 0, "100M", DEFAULT_NM, "", "");

        List<AlignData> alignments = Lists.newArrayList(alignment1, alignment2);

        List<AlignData> validAlignments = Lists.newArrayList();
        List<AlignData> lowQualAlignments = Lists.newArrayList();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(0, lowQualAlignments.size());
        assertEquals(2, validAlignments.size());
        assertEquals(100, validAlignments.get(0).adjustedAlignment());
        assertEquals(100, validAlignments.get(1).adjustedAlignment());
        assertEquals(60, validAlignments.get(0).modifiedMapQual(), 0.1);

        // now with lower alignment score and overlap
        AlignData zeroAlign = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 200), 0, 1,
                0, 1, 0, "1M", DEFAULT_NM, "", "");

        alignment1 = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 200), 0, 80,
                60, 80, 0, "110M", DEFAULT_NM, "", "");

        alignment2 = new AlignData(
                new ChrBaseRegion(CHR_2, 201, 300), 60, 140,
                100, 80, 0, "100M", DEFAULT_NM, "", "");

        AlignData alignment3 = new AlignData(
                new ChrBaseRegion(CHR_3, 1, 100), 130, 200,
                60, 70, 0, "100M", DEFAULT_NM, "", "");

        alignments = Lists.newArrayList(zeroAlign, alignment1, zeroAlign, alignment2, zeroAlign, alignment3);

        validAlignments.clear();
        lowQualAlignments.clear();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(4, lowQualAlignments.size());
        assertEquals(2, validAlignments.size());

        assertEquals(60, alignment1.adjustedAlignment());
        assertEquals(60, alignment3.adjustedAlignment());
        assertEquals(48, alignment2.adjustedAlignment());
    }

    @Test
    public void testAlignmentShortAdjustedAlignedLengthFilter()
    {
        AssemblyAlignment assemblyAlignment = createAssemblyAlignment(
                CHR_1, 200, FORWARD, CHR_2, 200, REVERSE, "");

        // test exclusion of short adjusted alignment-length segments
        AlignData alignment1 = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 200), 0, 80,
                60, 80, 0, "80M", DEFAULT_NM, "", "");

        String altAlignment = "2,+200000,100M,0";
        AlignData alignment2 = new AlignData(
                new ChrBaseRegion(CHR_2, 201, 220), 80, 100,
                60, 20, 0, "20M", DEFAULT_NM, altAlignment, "");

        AlignData alignment3 = new AlignData(
                new ChrBaseRegion(CHR_3, 1, 100), 100, 200,
                60, 100, 0, "100M", DEFAULT_NM, "", "");

        List<AlignData> validAlignments = Lists.newArrayList();
        List<AlignData> lowQualAlignments = Lists.newArrayList();

        List<AlignData> alignments = Lists.newArrayList(alignment1, alignment2, alignment3);

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(1, lowQualAlignments.size());
        assertEquals(2, validAlignments.size());

        // short alignments are permitted on the outside
        alignment1 = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 120), 0, 20,
                60, 80, 0, "20M", DEFAULT_NM, altAlignment, "");

        alignment2 = new AlignData(
                new ChrBaseRegion(CHR_2, 201, 300), 20, 100,
                60, 80, 0, "100M", DEFAULT_NM, "", "");

        alignment3 = new AlignData(
                new ChrBaseRegion(CHR_2, 350, 370), 100, 120,
                60, 100, 0, "20M", DEFAULT_NM, altAlignment, "");

        validAlignments.clear();
        lowQualAlignments.clear();

        alignments = Lists.newArrayList(alignment1, alignment2, alignment3);

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(3, validAlignments.size());
    }

    @Test
    public void testAlignmentLowMapQualFilters()
    {
        AssemblyAlignment assemblyAlignment = AssemblyTestUtils.createAssemblyAlignment(
                mRefGenome, CHR_1, 300, FORWARD, CHR_1, 350, REVERSE, "", "",
                50, 200);

        // scenario 1: based on chr7-comp 1c:
        String cigar = "100M";

        String altAlignment1 = "2,+200000,100M,0";

        AlignData alignment1 = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 200), 0, 100,
                0, 100, 0, cigar, DEFAULT_NM, altAlignment1, "");

        String altAlignment2 = "1,+2000,100M,0";

        AlignData alignment2 = new AlignData(
                new ChrBaseRegion(CHR_1, 20000, 20100), 101, 200,
                0, 100, 0, cigar, DEFAULT_NM, altAlignment2, "");

        List<AlignData> alignments = Lists.newArrayList(alignment1, alignment2);

        List<AlignData> validAlignments = Lists.newArrayList();
        List<AlignData> lowQualAlignments = Lists.newArrayList();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(0, lowQualAlignments.size());
        assertEquals(2, validAlignments.size());

        assertTrue(alignment1.hasLowMapQualAlignment());
        assertTrue(alignment2.hasLowMapQualAlignment());

        assertNull(alignment1.selectedAltAlignment());
        assertEquals(1, alignment1.unselectedAltAlignments().size());
        assertEquals(2000, alignment2.selectedAltAlignment().Position);
        assertTrue(alignment1.hasLowMapQualShortSvLink());
        assertTrue(alignment2.hasLowMapQualShortSvLink());

        // scenario 2: based on chr7-comp 1a:

        // 0 = {AlignData@3222} "7:125745443-125746123:1 681M1314S seq(0-681 adj=0-680) score(681) flags(0) mapQual(60 adj=50) aligned(681 adj=621)"
        //1 = {AlignData@3223} "7:126166901-126167444:1 544M seq(0-544 adj=682-1225) score(544) flags(0) mapQual(60 adj=52) aligned(544 adj=506)"
        //2 = {AlignData@3224} "7:143939546-143939814:-1 269M seq(0-269 adj=1231-1499) score(269) flags(16) mapQual(0 adj=0) aligned(269 adj=245)"
            // 7,+144005351,269M,0;
        //3 = {AlignData@3225} "7:143936039-143936532:-1 494M seq(0-494 adj=1501-1994) score(494) flags(16) mapQual(0 adj=0) aligned(494 adj=478)"
            // 7,+144008633,494M,0;

        alignment1 = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 200), 0, 100,
                60, 100, 0, cigar, DEFAULT_NM, altAlignment1, "");

        alignment2 = new AlignData(
                new ChrBaseRegion(CHR_1, 201, 300), 101, 200,
                60, 100, 0, cigar, DEFAULT_NM, altAlignment2, "");

        String altAlignment3 = "2,+200000,100M,0";

        AlignData alignment3 = new AlignData(
                new ChrBaseRegion(CHR_2, 101, 200), 201, 300,
                0, 100, 0, cigar, DEFAULT_NM, altAlignment3, "");

        String altAlignment4 = "1,+2000,100M,0";

        AlignData alignment4 = new AlignData(
                new ChrBaseRegion(CHR_2, 20000, 20100), 301, 400,
                0, 100, 0, cigar, DEFAULT_NM, altAlignment4, "");

        alignments = Lists.newArrayList(alignment1, alignment2, alignment3, alignment4);

        validAlignments.clear();
        lowQualAlignments.clear();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(0, lowQualAlignments.size());
        assertEquals(4, validAlignments.size());

        assertTrue(alignment3.hasLowMapQualAlignment());
        assertTrue(alignment4.hasLowMapQualAlignment());

        assertNull(alignment3.selectedAltAlignment());
        assertEquals(1, alignment3.unselectedAltAlignments().size());
        assertNull(alignment4.selectedAltAlignment());
        assertEquals(1, alignment4.unselectedAltAlignments().size());
        assertTrue(alignment3.hasLowMapQualShortSvLink());
        assertTrue(alignment4.hasLowMapQualShortSvLink());
   }

    @Test
    public void testKeepAdjacentLocalShortAdjustedAlignments()
    {
        AssemblyAlignment assemblyAlignment = AssemblyTestUtils.createAssemblyAlignment(
                mRefGenome, CHR_1, 300, FORWARD, CHR_1, 350, REVERSE, "", "",
                50, 200);

        // scenario 1: based on chr7-comp 1c:
        String cigar = "100M";
        String middleCigar = "50M";

        AlignData alignment1 = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 200), 0, 100,
                60, 100, 0, cigar, DEFAULT_NM, null, "");

        AlignData alignment2 = new AlignData(
                new ChrBaseRegion(CHR_1, 1101, 1150), 101, 150,
                60, 40, 0, middleCigar, DEFAULT_NM, null, "");

        AlignData alignment3 = new AlignData(
                new ChrBaseRegion(CHR_1, 20000, 20100), 201, 300,
                60, 100, 0, cigar, DEFAULT_NM, null, "");

        List<AlignData> alignments = Lists.newArrayList(alignment1, alignment2, alignment3);

        List<AlignData> validAlignments = Lists.newArrayList();
        List<AlignData> lowQualAlignments = Lists.newArrayList();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(0, lowQualAlignments.size());
        assertEquals(3, validAlignments.size());

        // rescued by upper alignment
        alignment2 = new AlignData(
                new ChrBaseRegion(CHR_1, 19101, 19150), 101, 150,
                60, 40, 0, middleCigar, DEFAULT_NM, null, "");

        alignments = Lists.newArrayList(alignment1, alignment2, alignment3);

        validAlignments.clear();
        lowQualAlignments.clear();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(0, lowQualAlignments.size());
        assertEquals(3, validAlignments.size());

        // no longer rescued if remote or too far away
        alignment2 = new AlignData(
                new ChrBaseRegion(CHR_2, 1101, 1150), 101, 150,
                60, 40, 0, middleCigar, DEFAULT_NM, null, "");

        alignments = Lists.newArrayList(alignment1, alignment2, alignment3);

        validAlignments.clear();
        lowQualAlignments.clear();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(1, lowQualAlignments.size());
        assertEquals(2, validAlignments.size());

        alignment2 = new AlignData(
                new ChrBaseRegion(CHR_1, 10000, 10049), 101, 150,
                60, 40, 0, middleCigar, DEFAULT_NM, null, "");

        alignments = Lists.newArrayList(alignment1, alignment2, alignment3);

        validAlignments.clear();
        lowQualAlignments.clear();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(1, lowQualAlignments.size());
        assertEquals(2, validAlignments.size());

        // only needs 1 anchoring high-qual alignment
        alignment2 = new AlignData(
                new ChrBaseRegion(CHR_1, 19101, 19150), 101, 150,
                60, 40, 0, middleCigar, DEFAULT_NM, null, "");

        alignments = Lists.newArrayList(alignment2, alignment3);

        validAlignments.clear();
        lowQualAlignments.clear();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(0, lowQualAlignments.size());
        assertEquals(2, validAlignments.size());

        // not rescued if has too many mismatches
        middleCigar = "20M10I20M";
        String mdTag = "10A10C10G10T10";

        alignment2 = new AlignData(
                new ChrBaseRegion(CHR_1, 19101, 19150), 101, 150,
                60, 40, 0, middleCigar, DEFAULT_NM, null, mdTag);

        alignments = Lists.newArrayList(alignment2, alignment3);

        validAlignments.clear();
        lowQualAlignments.clear();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(1, lowQualAlignments.size());
        assertEquals(1, validAlignments.size());
    }

    @Test
    public void testSpecificAlignmentHandling()
    {
        AssemblyAlignment assemblyAlignment = AssemblyTestUtils.createAssemblyAlignment(
                mRefGenome, CHR_1, 300, FORWARD, CHR_1, 350, REVERSE, "", "",
                50, 200);

        // scenario 1: based on chr7-comp 1c:
        String cigar = "100M";

        AlignData alignment1 = new AlignData(
                new ChrBaseRegion(CHR_1, 101, 200), 0, 100,
                60, 100, 0, cigar, DEFAULT_NM, "", "");

        int lowMapQual = SSX2_MAX_MAP_QUAL - 1;
        char altAlignmentOrientation = SSX2_GENE_ORIENT.asChar();

        ChrBaseRegion ssx2Region = SSX2_REGIONS_V37.get(0);

        String altAlignment2 = format("%s,%c%d,100M,%d",
                ssx2Region.Chromosome, altAlignmentOrientation, ssx2Region.start() + 1, lowMapQual);

        AlignData alignment2 = new AlignData(
                new ChrBaseRegion(CHR_2, 20000, 20100), 101, 200,
                0, 100, 0, cigar, DEFAULT_NM, altAlignment2, "");

        List<AlignData> alignments = Lists.newArrayList(alignment1, alignment2);

        List<AlignData> validAlignments = Lists.newArrayList();
        List<AlignData> lowQualAlignments = Lists.newArrayList();

        filterAlignments(assemblyAlignment, alignments, validAlignments, lowQualAlignments);

        assertEquals(0, lowQualAlignments.size());
        assertEquals(2, validAlignments.size());


        AlignData specificAlignment = validAlignments.stream().filter(x -> x.refLocation().overlaps(ssx2Region)).findFirst().orElse(null);
        assertNull(specificAlignment);

        assertTrue(alignment2.hasSelectedAltAlignment());

        AlternativeAlignment ssx2AltAlignment = alignment2.selectedAltAlignment();
        assertTrue(ssx2Region.containsPosition(ssx2AltAlignment.Chromosome, ssx2AltAlignment.Position));
    }

    private AssemblyAlignment createAssemblyAlignment(
            final String chrStart, int posStart, Orientation orientStart,
            final String chrEnd, int posEnd, Orientation orientEnd, final String insertedBases)
    {
        return AssemblyTestUtils.createAssemblyAlignment(
                mRefGenome, chrStart, posStart, orientStart, chrEnd, posEnd, orientEnd, insertedBases, "");
    }
}
