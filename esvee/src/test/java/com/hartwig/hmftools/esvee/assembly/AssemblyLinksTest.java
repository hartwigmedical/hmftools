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
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.getSupportTypeCount;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.isValidSupportCoordsVsJunction;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteRegionAssembler;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
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

        JunctionAssembly secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80);

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        posJunction = new Junction(CHR_1, 1000, FORWARD);
        negJunction = new Junction(CHR_1, 500, REVERSE);

        firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);
        secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80);

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

        Read juncRead2 = cloneRead(juncRead1, READ_ID_GENERATOR.nextId());

        JunctionAssembler junctionAssembler = new JunctionAssembler(firstJunction);
        JunctionAssembly firstAssembly = junctionAssembler.processJunction(List.of(juncRead1, juncRead2)).get(0);

        String juncReadBases2 = refSequence.substring(151, 201) + Nucleotides.reverseComplementBases(refSequence.substring(51, 101));

        Read juncRead3 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 151, juncReadBases2, "50M50S", CHR_1, 1, false);

        Read juncRead4 = cloneRead(juncRead3, READ_ID_GENERATOR.nextId());

        junctionAssembler = new JunctionAssembler(secondJunction);
        JunctionAssembly secondAssembly = junctionAssembler.processJunction(List.of(juncRead3, juncRead4)).get(0);

        Read discRead2 = createRead(
                juncRead2.id(), CHR_1, 120, refSequence.substring(120, 170), "50M", CHR_1, 51, false);
        discRead2.bamRecord().setSecondOfPairFlag(true);

        secondAssembly.addCandidateSupport(discRead2);

        // add the same junction read from the first assembly as a candidate for the second - this won't end up being counted as support
        // since it matches the exact read as support in the first assembly
        secondAssembly.addCandidateSupport(juncRead1);

        Read discRead1 = createRead(
                juncRead1.id(), CHR_1, 120, refSequence.substring(120, 170), "50M", CHR_1, 51, false);
        discRead1.bamRecord().setSecondOfPairFlag(true);
        secondAssembly.addCandidateSupport(discRead1);

        Read discRead3 = createRead(
                juncRead3.id(), CHR_1, 20, refSequence.substring(20, 70), "50M", CHR_1, 151, false);
        discRead3.bamRecord().setSecondOfPairFlag(true);

        firstAssembly.addCandidateSupport(discRead3);

        PhaseGroup phaseGroup = new PhaseGroup(firstAssembly, secondAssembly);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(refGenome, null, phaseGroup);

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
    public void testChainedLinks()
    {
        // similar to chromosome 3-10-12-3 in COLO - links are:
        // chr1:100:1 -> chr2:50:-1
        // chr2:50:-1 -> chr2:200:1
        // chr2:200:1 -> chr3:200:1
        // chr3:200:1 -> chr3:50:-1 (ie segment is reverse-linked)
        // chr3:50:-1 -> chr1:250:-1

        // NOTE: this is also the order of both the assemblies and links as they're added to the phase set


        String refSequence = REF_BASES_400;

        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refSequence);
        refGenome.RefGenomeMap.put(CHR_2, refSequence);
        refGenome.RefGenomeMap.put(CHR_3, refSequence);

        String refBases1 = REF_BASES_400.substring(1, 101);
        String extBases1 = refGenome.getBaseString(CHR_2, 50, 150);
        String assemblyBases1 = refBases1 + extBases1;

        String posJuncReadCigar = "100M50S";
        String negJuncReadCigar = "50S100M";
        String juncReadBases = REF_BASES_400.substring(0, 150); // unused but needs to be the correct length

        // chr1:100:1 - links to chr2:50:-1
        JunctionAssembly assembly1 = createAssembly(CHR_1, 100, FORWARD, assemblyBases1, refBases1.length() - 1);

        Read juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, juncReadBases, posJuncReadCigar, CHR_2, 100, true);

        assembly1.addJunctionRead(juncRead);

        // chr2:50:-1 - links to chr1:100:1 and faces chr2:200:1
        String refBases2 = REF_BASES_400.substring(50, 150);
        String extBases2 = refGenome.getBaseString(CHR_1, 51, 101);
        String assemblyBases2 = extBases2 + refBases2;

        JunctionAssembly assembly2 = createAssembly(CHR_2, 50, REVERSE, assemblyBases2, extBases2.length());

        assembly2.addJunctionRead(juncRead); // add prior so fragment IDs match & link

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 50, juncReadBases, negJuncReadCigar, CHR_3, 100, true);

        assembly2.addJunctionRead(juncRead);

        // chr2:200:1 - links to chr3:200:-1 inverted to faces chr2:50
        String refBases3 = REF_BASES_400.substring(101, 201);
        String extBases3 = refGenome.getBaseString(CHR_3, 151, 201);
        String assemblyBases3 = refBases3 + Nucleotides.reverseComplementBases(extBases3);

        JunctionAssembly assembly3 = createAssembly(CHR_2, 200, FORWARD, assemblyBases3, refBases3.length() - 1);

        assembly3.addJunctionRead(juncRead);

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 101, juncReadBases, posJuncReadCigar, CHR_3, 100, true);

        assembly3.addJunctionRead(juncRead);

        // chr3:200:1 - links to chr2:200:1 inverted and faces chr3:50:-1
        String refBases4 = REF_BASES_400.substring(101, 201);
        String extBases4 = refGenome.getBaseString(CHR_2, 151, 201);
        String assemblyBases4 = refBases4 + Nucleotides.reverseComplementBases(extBases4);

        JunctionAssembly assembly4 = createAssembly(CHR_3, 200, FORWARD, assemblyBases4, refBases4.length() - 1);

        assembly4.addJunctionRead(juncRead);
        assembly1.addJunctionRead(juncRead);

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_3, 101, juncReadBases, posJuncReadCigar, CHR_2, 100, true);

        assembly4.addJunctionRead(juncRead);

        // chr3:50:-1 - links to chr1:250:-1 inverted and faces chr3:200:1
        String refBases5 = REF_BASES_400.substring(50, 150);
        String extBases5 = refGenome.getBaseString(CHR_1, 250, 300);
        String assemblyBases5 = Nucleotides.reverseComplementBases(extBases5) + refBases5;

        JunctionAssembly assembly5 = createAssembly(CHR_3, 50, REVERSE, assemblyBases5, extBases5.length());

        assembly5.addJunctionRead(juncRead);

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_3, 50, juncReadBases, negJuncReadCigar, CHR_1, 300, true);

        assembly5.addJunctionRead(juncRead);

        // chr1:250:-1 - links to chr3:50:-1
        String refBases6 = REF_BASES_400.substring(250, 350);
        String extBases6 = refGenome.getBaseString(CHR_3, 50, 100);
        String assemblyBases6 = Nucleotides.reverseComplementBases(extBases6) + refBases6;

        JunctionAssembly assembly6 = createAssembly(CHR_1, 250, REVERSE, assemblyBases6, extBases6.length());

        assembly6.addJunctionRead(juncRead);

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 250, juncReadBases, negJuncReadCigar, CHR_3, 100, true);

        assembly6.addJunctionRead(juncRead);
        assembly4.addJunctionRead(juncRead); // for the facing link

        // add discordant reads as candidates to check they are assigned as support correctly
        String discCigar = "50M";
        Read discRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 20, REF_BASES_400.substring(20, 70), discCigar, CHR_3, 100, true);
        Read discRead2 = createRead(
                discRead1.id(), CHR_3, 100, REF_BASES_400.substring(100, 150), discCigar, CHR_1, 20, false);

        assertTrue(isValidSupportCoordsVsJunction(discRead1, assembly1.junction().isForward(), assembly1.junction().Position));
        assembly1.addCandidateSupport(discRead1);

        assertTrue(isValidSupportCoordsVsJunction(discRead2, assembly4.junction().isForward(), assembly4.junction().Position));
        assembly4.addCandidateSupport(discRead2);

        assertFalse(isValidSupportCoordsVsJunction(discRead2, assembly5.junction().isForward(), assembly5.junction().Position));
        assembly5.addCandidateSupport(discRead2);

        Read discRead3 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 100, REF_BASES_400.substring(100, 150), discCigar, CHR_1, 300, true);
        Read discRead4 = createRead(
                discRead3.id(), CHR_1, 300, REF_BASES_400.substring(300, 350), discCigar, CHR_2, 100, false);
        discRead4.bamRecord().setReadNegativeStrandFlag(true);

        assertFalse(isValidSupportCoordsVsJunction(discRead3, assembly2.junction().isForward(), assembly2.junction().Position));
        assertTrue(isValidSupportCoordsVsJunction(discRead3, assembly3.junction().isForward(), assembly3.junction().Position));
        assembly3.addCandidateSupport(discRead3);

        assertFalse(isValidSupportCoordsVsJunction(discRead4, assembly1.junction().isForward(), assembly1.junction().Position));
        assertTrue(isValidSupportCoordsVsJunction(discRead4, assembly6.junction().isForward(), assembly6.junction().Position));
        assembly6.addCandidateSupport(discRead4);

        PhaseGroup phaseGroup = new PhaseGroup(assembly1, assembly2);
        phaseGroup.addAssembly(assembly3);
        phaseGroup.addAssembly(assembly4);
        phaseGroup.addAssembly(assembly5);
        phaseGroup.addAssembly(assembly6);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(refGenome, new RemoteRegionAssembler(refGenome, null), phaseGroup);
        phaseSetBuilder.buildPhaseSets();

        assertEquals(1, phaseGroup.phaseSets().size());
        PhaseSet phaseSet = phaseGroup.phaseSets().get(0);

        assertEquals(5, phaseSet.assemblyLinks().size());

        hasAssemblyLink(phaseSet.assemblyLinks(), assembly1, assembly2, LinkType.SPLIT);
        hasAssemblyLink(phaseSet.assemblyLinks(), assembly2, assembly3, LinkType.FACING);
        hasAssemblyLink(phaseSet.assemblyLinks(), assembly3, assembly4, LinkType.SPLIT);
        hasAssemblyLink(phaseSet.assemblyLinks(), assembly4, assembly5, LinkType.FACING);
        hasAssemblyLink(phaseSet.assemblyLinks(), assembly5, assembly6, LinkType.SPLIT);

        assertEquals(3, assembly1.supportCount());
        assertEquals(2, assembly2.supportCount());
        assertEquals(3, assembly3.supportCount());
        assertEquals(4, assembly4.supportCount());
        assertEquals(2, assembly5.supportCount());
        assertEquals(3, assembly6.supportCount());
        assertEquals(1, getSupportTypeCount(assembly1, DISCORDANT));
        assertEquals(0, getSupportTypeCount(assembly2, DISCORDANT));
        assertEquals(1, getSupportTypeCount(assembly3, DISCORDANT));
        assertEquals(1, getSupportTypeCount(assembly4, DISCORDANT));
        assertEquals(0, getSupportTypeCount(assembly5, DISCORDANT));
        assertEquals(1, getSupportTypeCount(assembly6, DISCORDANT));
    }

    private static boolean hasAssemblyLink(
            final List<AssemblyLink> links, final JunctionAssembly assembly1, final JunctionAssembly assembly2, final LinkType linkType)
    {
        return links.stream().anyMatch(x -> x.type() == linkType && x.hasAssembly(assembly1) && x.hasAssembly(assembly2));
    }
}
