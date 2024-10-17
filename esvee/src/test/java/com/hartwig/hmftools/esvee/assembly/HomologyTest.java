package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.TestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.REF_GENOME;
import static com.hartwig.hmftools.esvee.alignment.HomologyData.determineHomology;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAlignment;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssemblyAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.esvee.alignment.AlignData;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.alignment.Breakend;
import com.hartwig.hmftools.esvee.alignment.BreakendBuilder;
import com.hartwig.hmftools.esvee.alignment.HomologyData;
import com.hartwig.hmftools.esvee.alignment.MdTag;
import com.hartwig.hmftools.esvee.alignment.MdTagType;
import com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

import org.junit.Test;

public class HomologyTest
{
    @Test
    public void testHomology()
    {
        AlignData alignmentStart = createAlignment(CHR_1, 100, 150, 0, 50, "51M");
        AlignData alignmentEnd = createAlignment(CHR_1, 100, 150, 51, 100, "50M");

        String assemblyOverlap = "";

        assertNull(determineHomology(assemblyOverlap, alignmentStart, alignmentEnd, REF_GENOME));

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

    @Test
    public void testPositionShift()
    {
        HomologyData homology = new HomologyData("", -2, 1, -2, 1);
        assertEquals(-1, homology.positionAdjustment(FORWARD));
        assertEquals(2, homology.positionAdjustment(REVERSE));
    }

    @Test
    public void testHomologyScenario()
    {
        MockRefGenome refGenome = new MockRefGenome();
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_400);

        // first a DEL with 3 bases of exact homology - 'TAG'
        String homologyBases = "TAG";
        AssemblyAlignment assemblyAlignment = createAssemblyAlignment(
                refGenome, CHR_1, 121, FORWARD, CHR_1, 202, REVERSE, "", homologyBases);

        BreakendBuilder breakendBuilder = new BreakendBuilder(refGenome, assemblyAlignment);

        AlignData alignment1 = createAlignment(CHR_1, 22, 121, 0, 100, "100M50S");
        AlignData alignment2 = createAlignment(CHR_1, 202, 302, 97, 197, "50S100M");

        List<AlignData> alignments = Lists.newArrayList(alignment1, alignment2);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        Breakend first = assemblyAlignment.breakends().get(0);
        assertEquals(DEL, first.svType());
        assertEquals(120, first.Position);
        assertEquals(FORWARD, first.Orient);
        assertNotNull(first.Homology);
        assertEquals(homologyBases, first.Homology.Homology);
        assertEquals(-2, first.Homology.ExactStart);
        assertEquals(1, first.Homology.ExactEnd);

        Breakend second = assemblyAlignment.breakends().get(1);
        assertEquals(204, second.Position);
        assertEquals(REVERSE, second.Orient);
        assertEquals(-2, second.Homology.ExactStart);
        assertEquals(1, second.Homology.ExactEnd);


        // test 2: positive INV
        assemblyAlignment = createAssemblyAlignment(
                refGenome, CHR_1, 121, FORWARD, CHR_1, 340, FORWARD, "", homologyBases);

        breakendBuilder = new BreakendBuilder(refGenome, assemblyAlignment);

        alignment1 = createAlignment(CHR_1, 22, 121, 0, 100, "100M50S");

        alignment2 = createAlignment(
                CHR_1, 241, 340, true, DEFAULT_MAP_QUAL, 100, 0, 100,
                "100M50S", "");

        alignments = Lists.newArrayList(alignment1, alignment2);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        first = assemblyAlignment.breakends().get(0);
        assertEquals(INV, first.svType());
        assertEquals(120, first.Position);
        assertEquals(FORWARD, first.Orient);
        assertNotNull(first.Homology);
        assertEquals(homologyBases, first.Homology.Homology);
        assertEquals(-2, first.Homology.ExactStart);
        assertEquals(1, first.Homology.ExactEnd);

        second = assemblyAlignment.breakends().get(1);
        assertEquals(338, second.Position);
        assertEquals(FORWARD, second.Orient);
        assertEquals("CTA", second.Homology.Homology);
        assertEquals(-1, second.Homology.ExactStart);
        assertEquals(2, second.Homology.ExactEnd);


        // test 3: negative INV
        int pos1 = 176;
        int pos2 = 280;

        Junction junction1 = new Junction(CHR_1, pos1, REVERSE);
        Junction junction2 = new Junction(CHR_1, pos2, REVERSE);

        int refLength = 100;
        int extensionLength = 50;

        String firstRefBases = REF_BASES_400.substring(pos1, pos1 + refLength); // starts with TAG
        String secondRefBases = REF_BASES_400.substring(pos2, pos2 + refLength); // runs into CTA

        String firstExtensionBases = Nucleotides.reverseComplementBases(secondRefBases.substring(3, 3 + extensionLength)); // first 80 bases of second's ref, exact match and no insert
        String firstAssemblyBases = firstExtensionBases + firstRefBases;

        byte[] baseQuals = buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(junction1, firstAssemblyBases.getBytes(), baseQuals, extensionLength);

        String secondExtensionBases =  Nucleotides.reverseComplementBases(firstRefBases.substring(3, 3 + extensionLength)); // last 80 bases of first's ref, again an exact match
        String secondAssemblyBases = secondExtensionBases + secondRefBases;

        JunctionAssembly secondAssembly = new JunctionAssembly(junction2, secondAssemblyBases.getBytes(), baseQuals, extensionLength);

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        // order passed in doesn't matter
        homologyBases = "CTA"; // for the second breakend, which means last 3 ref bases of first breakend reads TAG
        String homologyBasesRev = Nucleotides.reverseComplementBases(homologyBases);
        AssemblyLink link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);
        assertNotNull(link);
        assertEquals(homologyBasesRev, link.overlapBases());

        assemblyAlignment = createAssemblyAlignment(link);

        breakendBuilder = new BreakendBuilder(refGenome, assemblyAlignment);

        alignment1 = createAlignment(
                CHR_1, pos1, pos1 + refLength - 1, true, DEFAULT_MAP_QUAL, 100, 97, 197,
                "50S100M", "");

        alignment2 = createAlignment(CHR_1, pos2, pos2 + refLength - 1, 97, 197, "50S100M");

        alignments = Lists.newArrayList(alignment1, alignment2);

        breakendBuilder.formBreakends(alignments);

        assertEquals(2, assemblyAlignment.breakends().size());

        first = assemblyAlignment.breakends().get(0);
        assertEquals(INV, first.svType());
        assertEquals(pos1 + 2, first.Position);
        assertEquals(REVERSE, first.Orient);
        assertNotNull(first.Homology);
        assertEquals(homologyBasesRev, first.Homology.Homology);
        assertEquals(-2, first.Homology.ExactStart);
        assertEquals(1, first.Homology.ExactEnd);

        second = assemblyAlignment.breakends().get(1);
        assertEquals(pos2 + 1, second.Position);
        assertEquals(REVERSE, second.Orient);
        assertEquals(homologyBases, second.Homology.Homology);
        assertEquals(-1, second.Homology.ExactStart);
        assertEquals(2, second.Homology.ExactEnd);
    }

    @Test
    public void testMdTagParsing()
    {
        MdTag mdTag = new MdTag("5G3");

        assertEquals(3, mdTag.elements().size());
        assertEquals(9, mdTag.baseLength());
        int index = 0;
        assertEquals(5, mdTag.elements().get(index).Length);
        assertEquals(MdTagType.MATCH, mdTag.elements().get(index).Type);

        assertEquals(MdTagType.SNV, mdTag.elements().get(++index).Type);
        assertEquals(3, mdTag.elements().get(++index).Length);
        assertEquals(MdTagType.MATCH, mdTag.elements().get(index).Type);

        index = 0;
        mdTag = new MdTag("0A0A3T");
        assertEquals(4, mdTag.elements().size());
        assertEquals(6, mdTag.baseLength());
        assertEquals(MdTagType.SNV, mdTag.elements().get(index).Type);

        assertEquals(MdTagType.SNV, mdTag.elements().get(++index).Type);
        assertEquals(3, mdTag.elements().get(++index).Length);
        assertEquals(MdTagType.SNV, mdTag.elements().get(++index).Type);

        // deletes
        index = 0;
        mdTag = new MdTag("5^AA0T5^GGG5");
        assertEquals(6, mdTag.elements().size());
        assertEquals(21, mdTag.baseLength());
        assertEquals(MdTagType.MATCH, mdTag.elements().get(index).Type);
        assertEquals(5, mdTag.elements().get(index).Length);

        assertEquals(MdTagType.DEL, mdTag.elements().get(++index).Type);
        assertEquals(2, mdTag.elements().get(index).Length);

        assertEquals(MdTagType.SNV, mdTag.elements().get(++index).Type);

        assertEquals(MdTagType.MATCH, mdTag.elements().get(++index).Type);
        assertEquals(5, mdTag.elements().get(index).Length);

        assertEquals(MdTagType.DEL, mdTag.elements().get(++index).Type);
        assertEquals(3, mdTag.elements().get(index).Length);

        assertEquals(MdTagType.MATCH, mdTag.elements().get(++index).Type);
        assertEquals(5, mdTag.elements().get(index).Length);

        String mdTagStr = "12G16T12G1T8C0A6C1C0C9C0C0T0T7T20A0C0C0C7C0A1C5^G1T8A9G3C9C0T8C3A6T3^CTC14C9G13C679";

        mdTag = new MdTag(mdTagStr);

        assertEquals(61, mdTag.elements().size());
    }

    @Test
    public void testMdSequenceComparisons()
    {
        // String seq1 = "AAAAAAAA";
        // String seq2 = "AAAACAAA";

        /*
        String mdTag1 = "8";
        String mdTag2 = "4C3";

        // byte[] mdSeq1 = MdTagComparer.extractSequence(mdTag1, seq1, 0, seq1.length() - 1);
        //byte[] mdSeq2 = MdTagComparer.extractSequence(mdTag2, seq2, 0, seq1.length() - 1);

        int endMismatchIndex = compareMdTagSequences(
                mdTag1, 0, 7, mdTag2, 0, 7);

        assertEquals(3, endMismatchIndex);

        mdTag1 = "12G16T12G1T8C0A6C1C0C9C0C0T0T7T20A0C0C0C7C0A1C5^G1T8A9G3C9C0T8C3A6T3^CTC14C9G13C679";
        mdTag2 = "76C644";

        int comparisonLength = 244;
        int secondLength = 721;

        endMismatchIndex = compareMdTagSequences(
                mdTag1, 0, comparisonLength - 1, mdTag2,
                secondLength - comparisonLength, secondLength - 1);

        assertEquals(3, endMismatchIndex);
        */


    }
}
