package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_600;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.getSupportTypeCount;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteReadExtractor;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;

import org.junit.Test;

public class AssemblyLinksTest
{
    private static final String INSERTED_BASES = "GGGGGGG";

    @Test
    public void testAssemblyBasicSplits()
    {
        String firstRefBases = REF_BASES_200.substring(0, 100);
        String secondRefBases = REF_BASES_200.substring(100, 200);

        // first a basic exact match junction
        Junction posJunction = new Junction(CHR_1, 100, FORWARD);
        Junction negJunction = new Junction(CHR_2, 200, REVERSE);

        int extensionLength = 80;
        String firstExtensionBases =  secondRefBases.substring(0, extensionLength); // first 80 bases of second's ref, exact match and no insert
        String firstAssemblyBases = firstRefBases + firstExtensionBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        String secondExtensionBases =  firstRefBases.substring(20, 100); // last 80 bases of first's ref, again an exact match
        String secondAssemblyBases = secondExtensionBases + secondRefBases;

        JunctionAssembly secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80);

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        // order passed in doesn't matter
        AssemblyLink link = assemblyLinker.tryAssemblyOverlap(secondAssembly, firstAssembly);
        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(BND, link.svType());


        // test 2: DEL with inserted bases
        posJunction = new Junction(CHR_1, 500, FORWARD);
        negJunction = new Junction(CHR_1, 1000, REVERSE);

        int insertedBaseLength = INSERTED_BASES.length();
        firstAssemblyBases = firstRefBases + INSERTED_BASES + firstExtensionBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        secondAssemblyBases = secondExtensionBases + INSERTED_BASES + secondRefBases;

        secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80 + insertedBaseLength);

        link = assemblyLinker.tryAssemblyOverlap(secondAssembly, firstAssembly);
        assertNotNull(link);
        assertEquals(INSERTED_BASES, link.insertedBases());
        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(DEL, link.svType());
    }

    @Test
    public void testAssemblyDups()
    {
        String firstRefBases = REF_BASES_200.substring(0, 100);
        String secondRefBases = REF_BASES_200.substring(100, 200);

        // first a basic exact match junction
        Junction posJunction = new Junction(CHR_1, 100, FORWARD);
        Junction negJunction = new Junction(CHR_2, 200, REVERSE);

        int extensionLength = 80;
        String firstExtensionBases =
                secondRefBases.substring(0, extensionLength); // first 80 bases of second's ref, exact match and no insert
        String firstAssemblyBases = firstRefBases + firstExtensionBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        String secondExtensionBases = firstRefBases.substring(20, 100); // last 80 bases of first's ref, again an exact match
        String secondAssemblyBases = secondExtensionBases + secondRefBases;

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        posJunction = new Junction(CHR_1, 1000, FORWARD);
        negJunction = new Junction(CHR_1, 500, REVERSE);

        firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);
        JunctionAssembly secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80);

        AssemblyLink link = assemblyLinker.tryAssemblyOverlap(secondAssembly, firstAssembly);
        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(DUP, link.svType());

        // now with a mismatch in the overlapping bases in the second assembly - has to be close to the junction point
        String mismatchSecondExtensionBases = secondExtensionBases.substring(0, 60)
                + MockRefGenome.getNextBase(secondExtensionBases.substring(60, 61)) + secondExtensionBases.substring(61, extensionLength);

        secondAssemblyBases = mismatchSecondExtensionBases + secondRefBases;
        secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80);

        link = assemblyLinker.tryAssemblyOverlap(secondAssembly, firstAssembly);
        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(DUP, link.svType());

        // DUP with overlapping bases - in this case the +ve assembly's extension bases start 20 bases into the ref sequence of the -ve
        // first assembly: ref bases + shared overlap + extension (being 20 bases into the ref bases of the second)
        // second assembly: extension (being ref bases 20-80 of first) + overlap + ref bases

        posJunction = new Junction(CHR_1, 1000, FORWARD);
        negJunction = new Junction(CHR_1, 500, REVERSE);

        extensionLength = 60;

        int overlapLength = 20;
        String overlapBases = secondRefBases.substring(0, overlapLength); // first 20 ref bases are shared
        firstRefBases = REF_BASES_200.substring(0, 80) + overlapBases;
        firstExtensionBases = secondRefBases.substring(overlapLength, overlapLength + extensionLength);
        firstAssemblyBases = firstRefBases + firstExtensionBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        secondExtensionBases = firstRefBases.substring(overlapLength, overlapLength + extensionLength);
        secondAssemblyBases = secondExtensionBases + secondRefBases;

        secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 60);

        link = assemblyLinker.tryAssemblyOverlap(secondAssembly, firstAssembly);
        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(DUP, link.svType());
        assertEquals(overlapBases, link.overlapBases());
    }

    @Test
    public void testAssemblyPositiveInversionSplits()
    {
        String firstRefBases = REF_BASES_200.substring(0, 100);
        String secondRefBases = REF_BASES_200.substring(100, 200);

        // local inversion
        Junction firstJunction = new Junction(CHR_1, 100, FORWARD);
        Junction secondJunction = new Junction(CHR_1, 200, FORWARD);

        String firstExtensionBases = reverseComplementBases(secondRefBases.substring(20, 100));
        String firstAssemblyBases = firstRefBases + firstExtensionBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        String secondExtensionBases = reverseComplementBases(firstRefBases.substring(20, 100));
        String secondAssemblyBases = secondRefBases + secondExtensionBases;

        JunctionAssembly secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, 99);

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        AssemblyLink link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);

        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());

        // test 2: same again but with an insert sequence
        firstAssemblyBases = firstRefBases + INSERTED_BASES + firstExtensionBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        secondAssemblyBases = secondRefBases + reverseComplementBases(INSERTED_BASES) + secondExtensionBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, 99);

        link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);
        assertNotNull(link);
        assertEquals(INSERTED_BASES, link.insertedBases());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());


        // test 3: 1 base of homology
        //                          10        20        30        40        50        60        70        80        90        100
        //                01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        firstRefBases =  "CAAACAAAAATTAAAACAAAAAGGCTAAATACTGCGTGATTTTGTTTACACGGAATTCTAAAAAGAAAAAAAAAAAAAAGACAACACAGAGAAAAAAACA";
        secondRefBases = "AATTCAAAGGACTGGCTTTAGCCTCCAGCTACTGAACAGCTTTGGTTAATAAGAAGCACCTGCAGAAGACTAGAAGACAACAGGAAAGGGGTTGTCTAGT";

        firstExtensionBases = reverseComplementBases(secondRefBases.substring(19, 99));
        firstAssemblyBases = firstRefBases + firstExtensionBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        secondExtensionBases = reverseComplementBases(firstRefBases.substring(19, 99));
        secondAssemblyBases = secondRefBases + secondExtensionBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, 99);

        link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);

        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());
        assertEquals("A", link.overlapBases());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());

        // test 4: again but with a base from the second's extension bases, so the mismatch must be factored into the offset

        secondExtensionBases = reverseComplementBases(firstRefBases.substring(19, 62) + firstRefBases.substring(63, 99));
        secondAssemblyBases = secondRefBases + secondExtensionBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, 99);

        link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);

        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());
        assertEquals("A", link.overlapBases());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());
    }

    @Test
    public void testAssemblyNegativeInvertedSplits()
    {
        String firstRefBases = REF_BASES_400.substring(0, 100);
        String secondRefBases = REF_BASES_400.substring(100, 200);

        // a negative inversion
        Junction firstJunction = new Junction(CHR_1, 100, REVERSE);
        Junction secondJunction = new Junction(CHR_1, 200, REVERSE);

        int extensionLength = 80;
        String firstExtensionBases = reverseComplementBases(secondRefBases.substring(0, extensionLength)); // the first 80 bases
        String firstAssemblyBases = firstExtensionBases + firstRefBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, extensionLength);

        String secondExtensionBases = reverseComplementBases(firstRefBases.substring(0, extensionLength));
        String secondAssemblyBases = secondExtensionBases + secondRefBases;

        JunctionAssembly secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, extensionLength);

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        AssemblyLink link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);
        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());

        // test 2: same again but with an insert sequence - always taken from the first assembly without reverse-complimenting
        int insertedBaseLength = INSERTED_BASES.length();
        firstAssemblyBases = firstExtensionBases + INSERTED_BASES + firstRefBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, extensionLength + insertedBaseLength);

        secondAssemblyBases = secondExtensionBases + reverseComplementBases(INSERTED_BASES) + secondRefBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, extensionLength + insertedBaseLength);

        link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);
        assertNotNull(link);
        assertEquals(INSERTED_BASES, link.insertedBases());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());

        // test 3: with a 5-base overlap
        int overlapLength = 5;
        int newExtensionLength = extensionLength - overlapLength;
        firstAssemblyBases = firstExtensionBases.substring(0, extensionLength - overlapLength) + firstRefBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, newExtensionLength);

        secondAssemblyBases = secondExtensionBases.substring(0, extensionLength - overlapLength) + secondRefBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, newExtensionLength);

        link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);
        assertNotNull(link);

        String overlapBases = firstRefBases.substring(0, overlapLength);
        assertEquals(overlapBases, link.overlapBases());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());
    }

    @Test
    public void testInversionSupport()
    {
        String refSequence = REF_BASES_400;

        MockRefGenome refGenome = new MockRefGenome(false);
        refGenome.RefGenomeMap.put(CHR_1, refSequence);

        // local inversion
        Junction firstJunction = new Junction(CHR_1, 100, FORWARD);
        Junction secondJunction = new Junction(CHR_1, 200, FORWARD);

        // create a read which is a discordant candidate for one junction and a junction read for the other
        String juncReadBases = refSequence.substring(51, 101) + Nucleotides.reverseComplementBases(refSequence.substring(151, 201));

        Read juncRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 51, juncReadBases, "50M50S", CHR_1, 1, false);

        Read juncRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 52, juncReadBases.substring(1), "49M50S", CHR_1, 1, false);

        JunctionAssembler junctionAssembler = new JunctionAssembler(firstJunction);
        JunctionAssembly firstAssembly = junctionAssembler.processJunction(List.of(juncRead1, juncRead2)).get(0);

        String juncReadBases2 = refSequence.substring(151, 201) + Nucleotides.reverseComplementBases(refSequence.substring(51, 101));

        Read juncRead3 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 151, juncReadBases2, "50M50S", CHR_1, 1, false);

        Read juncRead4 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 152, juncReadBases2.substring(1), "49M50S", CHR_1, 1, false);

        junctionAssembler = new JunctionAssembler(secondJunction);
        JunctionAssembly secondAssembly = junctionAssembler.processJunction(List.of(juncRead3, juncRead4)).get(0);

        Read discRead2 = createRead(
                juncRead2.id(), CHR_1, 120, refSequence.substring(120, 175), "55M", CHR_1, 51, false);
        discRead2.bamRecord().setSecondOfPairFlag(true);

        secondAssembly.addCandidateSupport(discRead2);

        // add the same junction read from the first assembly as a candidate for the second - this won't end up being counted as support
        // since it matches the exact read as support in the first assembly
        secondAssembly.addCandidateSupport(juncRead1);

        Read discRead1 = createRead(
                juncRead1.id(), CHR_1, 120, refSequence.substring(120, 175), "55M", CHR_1, 51, false);
        discRead1.bamRecord().setSecondOfPairFlag(true);
        secondAssembly.addCandidateSupport(discRead1);

        Read discRead3 = createRead(
                juncRead3.id(), CHR_1, 20, refSequence.substring(20, 80), "60M", CHR_1, 151, false);
        discRead3.bamRecord().setSecondOfPairFlag(true);

        firstAssembly.addCandidateSupport(discRead3);

        PhaseGroup phaseGroup = new PhaseGroup(firstAssembly, secondAssembly);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(
                refGenome, new RemoteReadExtractor(null), phaseGroup);

        phaseSetBuilder.buildPhaseSets();

        assertEquals(1, phaseGroup.phaseSets().size());
        PhaseSet phaseSet = phaseGroup.phaseSets().get(0);

        assertEquals(2, phaseSet.assemblies().size());
        AssemblyLink link = phaseSet.assemblyLinks().get(0);

        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());

        assertEquals(3, firstAssembly.supportCount());
        assertEquals(2, getSupportTypeCount(firstAssembly, JUNCTION));
        assertEquals(1, getSupportTypeCount(firstAssembly, DISCORDANT));

        assertEquals(4, secondAssembly.supportCount());
        assertEquals(2, getSupportTypeCount(secondAssembly, JUNCTION));
        assertEquals(2, getSupportTypeCount(secondAssembly, DISCORDANT));
    }

    @Test
    public void testExtBaseOnlyMatchedLinks()
    {
        // two assemblies only match in their extension bases, so need to form the full insert sequence drawing from both

        String assemblyRefBases = REF_BASES_400.substring(1, 101);
        String assemblyExtBases = REF_BASES_200.substring(0, 100);
        String assemblyBases = assemblyRefBases + assemblyExtBases;

        JunctionAssembly assembly1 = createAssembly(CHR_1, 100, FORWARD, assemblyBases, assemblyRefBases.length() - 1);

        String assemblyRefBases2 = REF_BASES_400.substring(300, 400);
        String assemblyExtBases2 = REF_BASES_200.substring(30, 130);
        String assemblyBases2 = assemblyExtBases2 + assemblyRefBases2;

        JunctionAssembly assembly2 = createAssembly(CHR_2, 300, REVERSE, assemblyBases2, assemblyExtBases2.length());

        AssemblyLink link = AssemblyLinker.tryAssemblyOverlap(assembly1, assembly2);
        assertNotNull(link);
        assertEquals(130, link.insertedBases().length());
        assertEquals(REF_BASES_200.substring(0, 130), link.insertedBases());

        // cannot form a link purely in extension bases if the assemblies form a local DEL or DUP
        assembly2 = createAssembly(CHR_1, 300, REVERSE, assemblyBases2, assemblyExtBases2.length());

        link = AssemblyLinker.tryAssemblyOverlap(assembly1, assembly2);
        assertNull(link);
    }

    @Test
    public void testFacingLinkConditions()
    {
        String assemblyRefBases = REF_BASES_400.substring(151, 251);
        String assemblyExtBases = REF_BASES_600.substring(1, 50);
        String assemblyBases = assemblyRefBases + assemblyExtBases;

        JunctionAssembly assemblyPos = createAssembly(CHR_1, 250, FORWARD, assemblyBases, assemblyRefBases.length() - 1);

        String assemblyRefBases2 = REF_BASES_400.substring(100, 200);
        String assemblyExtBases2 = REF_BASES_600.substring(1, 50);
        String assemblyBases2 = assemblyExtBases2 + assemblyRefBases2;

        JunctionAssembly assemblyNeg = createAssembly(CHR_1, 100, REVERSE, assemblyBases2, assemblyExtBases2.length());

        // no shared read
        AssemblyLink link = AssemblyLinker.tryAssemblyFacing(assemblyNeg, assemblyPos, Collections.emptyList());
        assertNull(link);

        Read juncRead1a = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 151, assemblyBases, "100M50S", CHR_1, 100, false);
        juncRead1a.bamRecord().setReadNegativeStrandFlag(true);

        assemblyPos.addJunctionRead(juncRead1a);

        Read juncRead1b = createRead(
                juncRead1a.id(), CHR_1, 100, assemblyBases2, "50S100M", CHR_1, 151, true);

        assemblyNeg.addJunctionRead(juncRead1b);

        link = AssemblyLinker.tryAssemblyFacing(assemblyNeg, assemblyPos, Collections.emptyList());
        assertNotNull(link);

        // invalid if too close or far away
        assemblyNeg = createAssembly(CHR_1, 225, REVERSE, assemblyBases2, assemblyExtBases2.length());

        juncRead1b = createRead(
                juncRead1a.id(), CHR_1, 225, assemblyBases2, "50S100M", CHR_1, 151, true);

        assemblyNeg.addJunctionRead(juncRead1b);

        link = AssemblyLinker.tryAssemblyFacing(assemblyNeg, assemblyPos, Collections.emptyList());
        assertNull(link);

        // different chromosome
        assemblyNeg = createAssembly(CHR_2, 100, REVERSE, assemblyBases2, assemblyExtBases2.length());

        juncRead1b = createRead(
                juncRead1a.id(), CHR_2, 100, assemblyBases2, "50S100M", CHR_1, 151, true);

        assemblyNeg.addJunctionRead(juncRead1b);

        link = AssemblyLinker.tryAssemblyFacing(assemblyNeg, assemblyPos, Collections.emptyList());
        assertNull(link);

        // matched reads have incorrect orientation
        assemblyNeg = createAssembly(CHR_1, 100, REVERSE, assemblyBases2, assemblyExtBases2.length());

        juncRead1b = createRead(
                juncRead1a.id(), CHR_1, 100, assemblyBases2, "50S100M", CHR_1, 151, true);

        assemblyNeg.addJunctionRead(juncRead1b);

        link = AssemblyLinker.tryAssemblyFacing(assemblyNeg, assemblyPos, Collections.emptyList());
        assertNotNull(link);

        assemblyNeg.support().clear();

        juncRead1b.bamRecord().setReadNegativeStrandFlag(true);
        assemblyNeg.addJunctionRead(juncRead1b);

        link = AssemblyLinker.tryAssemblyFacing(assemblyNeg, assemblyPos, Collections.emptyList());
        assertNull(link);

        // or reads facing away
        assemblyPos.support().clear();

        juncRead1a.bamRecord().setReadNegativeStrandFlag(false);
        assemblyPos.addJunctionRead(juncRead1a);

        link = AssemblyLinker.tryAssemblyFacing(assemblyNeg, assemblyPos, Collections.emptyList());
        assertNull(link);
    }
}
