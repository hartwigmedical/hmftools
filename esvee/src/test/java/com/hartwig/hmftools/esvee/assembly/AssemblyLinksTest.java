package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.createAssembly;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.formTestRefSequence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder;
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
        Junction posJunction = new Junction(CHR_1, 100, POS_ORIENT);
        Junction negJunction = new Junction(CHR_2, 200, NEG_ORIENT);

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

        // DEL with inserted bases
        posJunction = new Junction(CHR_1, 500, POS_ORIENT);
        negJunction = new Junction(CHR_1, 1000, NEG_ORIENT);

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
        Junction posJunction = new Junction(CHR_1, 100, POS_ORIENT);
        Junction negJunction = new Junction(CHR_2, 200, NEG_ORIENT);

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

        posJunction = new Junction(CHR_1, 1000, POS_ORIENT);
        negJunction = new Junction(CHR_1, 500, NEG_ORIENT);

        firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);
        secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80);

        // order passed in doesn't matter
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

        posJunction = new Junction(CHR_1, 1000, POS_ORIENT);
        negJunction = new Junction(CHR_1, 500, NEG_ORIENT);

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

        // order passed in doesn't matter
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
        Junction firstJunction = new Junction(CHR_1, 100, POS_ORIENT);
        Junction secondJunction = new Junction(CHR_1, 200, POS_ORIENT);

        String firstExtensionBases = reverseComplementBases(secondRefBases.substring(20, 100));
        String firstAssemblyBases = firstRefBases + firstExtensionBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        String secondExtensionBases = reverseComplementBases(firstRefBases.substring(20, 100));
        String secondAssemblyBases = secondRefBases + secondExtensionBases;

        JunctionAssembly secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, 99);

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        // order passed in doesn't matter
        AssemblyLink link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);
        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());

        // same again but with an insert sequence
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
    }

    @Test
    public void testAssemblyNegativeInvertedSplits()
    {
        String firstRefBases = REF_BASES_200.substring(0, 100);
        String secondRefBases = REF_BASES_200.substring(100, 200);

        // a negative inversion
        Junction firstJunction = new Junction(CHR_1, 100, NEG_ORIENT);
        Junction secondJunction = new Junction(CHR_1, 200, NEG_ORIENT);

        int extensionLength = 80;
        String firstExtensionBases = reverseComplementBases(secondRefBases.substring(0, extensionLength)); // the first 80 bases
        String firstAssemblyBases = firstExtensionBases + firstRefBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, extensionLength);

        String secondExtensionBases = reverseComplementBases(firstRefBases.substring(0, extensionLength));
        String secondAssemblyBases = secondExtensionBases + secondRefBases;

        JunctionAssembly secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, extensionLength);

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        // order passed in doesn't matter
        AssemblyLink link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);
        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());

        // same again but with an insert sequence - is it arbitrary who the bases appear reversed for?
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
        JunctionAssembly assembly1 = createAssembly(CHR_1, 100, POS_ORIENT, assemblyBases1, refBases1.length() - 1);

        Read juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, juncReadBases, posJuncReadCigar, CHR_2, 100, true);

        assembly1.addJunctionRead(juncRead);

        // chr2:50:-1 - links to chr1:100:1 and faces chr2:200:1
        String refBases2 = REF_BASES_400.substring(50, 150);
        String extBases2 = refGenome.getBaseString(CHR_1, 51, 101);
        String assemblyBases2 = extBases2 + refBases2;

        JunctionAssembly assembly2 = createAssembly(CHR_2, 50, NEG_ORIENT, assemblyBases2, extBases2.length());

        assembly2.addJunctionRead(juncRead); // add prior so fragment IDs match & link

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 50, juncReadBases, negJuncReadCigar, CHR_3, 100, true);

        assembly2.addJunctionRead(juncRead);

        // chr2:200:1 - links to chr3:200:-1 inverted to faces chr2:50
        String refBases3 = REF_BASES_400.substring(101, 201);
        String extBases3 = refGenome.getBaseString(CHR_3, 151, 201);
        String assemblyBases3 = refBases3 + Nucleotides.reverseComplementBases(extBases3);

        JunctionAssembly assembly3 = createAssembly(CHR_2, 200, POS_ORIENT, assemblyBases3, refBases3.length() - 1);

        assembly3.addJunctionRead(juncRead);

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 101, juncReadBases, posJuncReadCigar, CHR_3, 100, true);

        assembly3.addJunctionRead(juncRead);

        // chr3:200:1 - links to chr2:200:1 inverted and faces chr3:50:-1
        String refBases4 = REF_BASES_400.substring(101, 201);
        String extBases4 = refGenome.getBaseString(CHR_2, 151, 201);
        String assemblyBases4 = refBases4 + Nucleotides.reverseComplementBases(extBases4);

        JunctionAssembly assembly4 = createAssembly(CHR_3, 200, POS_ORIENT, assemblyBases4, refBases4.length() - 1);

        assembly4.addJunctionRead(juncRead);
        assembly1.addJunctionRead(juncRead);

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_3, 101, juncReadBases, posJuncReadCigar, CHR_2, 100, true);

        assembly4.addJunctionRead(juncRead);

        // chr3:50:-1 - links to chr1:250:-1 inverted and faces chr3:200:1
        String refBases5 = REF_BASES_400.substring(50, 150);
        String extBases5 = refGenome.getBaseString(CHR_1, 250, 300);
        String assemblyBases5 = Nucleotides.reverseComplementBases(extBases5) + refBases5;

        JunctionAssembly assembly5 = createAssembly(CHR_3, 50, NEG_ORIENT, assemblyBases5, extBases5.length());

        assembly5.addJunctionRead(juncRead);

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_3, 50, juncReadBases, negJuncReadCigar, CHR_1, 300, true);

        assembly5.addJunctionRead(juncRead);

        // chr1:250:-1 - links to chr3:50:-1
        String refBases6 = REF_BASES_400.substring(250, 350);
        String extBases6 = refGenome.getBaseString(CHR_3, 50, 100);
        String assemblyBases6 = Nucleotides.reverseComplementBases(extBases6) + refBases6;

        JunctionAssembly assembly6 = createAssembly(CHR_1, 250, NEG_ORIENT, assemblyBases6, extBases6.length());

        assembly6.addJunctionRead(juncRead);

        juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 250, juncReadBases, negJuncReadCigar, CHR_3, 100, true);

        assembly6.addJunctionRead(juncRead);
        assembly4.addJunctionRead(juncRead); // for the facing link

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
    }

    private static boolean hasAssemblyLink(
            final List<AssemblyLink> links, final JunctionAssembly assembly1, final JunctionAssembly assembly2, final LinkType linkType)
    {
        return links.stream().anyMatch(x -> x.type() == linkType && x.hasAssembly(assembly1) && x.hasAssembly(assembly2));
    }
}
