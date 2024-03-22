package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.types.AssemblyLink;
import com.hartwig.hmftools.esvee.types.Junction;
import com.hartwig.hmftools.esvee.types.JunctionAssembly;

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

        assertEquals(79, link.firstJunctionIndexInSecond());
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
        assertEquals(79, link.firstJunctionIndexInSecond());
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

        assertEquals(79, link.firstJunctionIndexInSecond());
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

        assertEquals(79, link.firstJunctionIndexInSecond());
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

        assertEquals(79, link.firstJunctionIndexInSecond());
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

        assertEquals(100, link.firstJunctionIndexInSecond());
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

        assertEquals(107, link.firstJunctionIndexInSecond());
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

        assertEquals(79, link.firstJunctionIndexInSecond());
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

        assertEquals(79, link.firstJunctionIndexInSecond());
        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());
    }
}
