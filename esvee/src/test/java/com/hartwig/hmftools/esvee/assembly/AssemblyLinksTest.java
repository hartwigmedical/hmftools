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

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.common.AssemblyLink;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;

import org.junit.Test;

public class AssemblyLinksTest
{
    @Test
    public void testAssemblyBasicSplits()
    {
        String firstRefBases =  REF_BASES_200.substring(0, 100);
        String secondRefBases =  REF_BASES_200.substring(100, 200);

        // first a basic exact match junction
        Junction posJunction = new Junction(CHR_1, 100, POS_ORIENT);
        Junction negJunction = new Junction(CHR_2, 200, NEG_ORIENT);

        String firstExtensionBases =  secondRefBases.substring(0, 80); // first 80 bases of second's ref, exact match and no insert
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

        // similar case but as a DUP
        posJunction = new Junction(CHR_1, 1000, POS_ORIENT);
        negJunction = new Junction(CHR_1, 500, NEG_ORIENT);

        firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);
        secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80);

        // order passed in doesn't matter
        link = assemblyLinker.tryAssemblyOverlap(secondAssembly, firstAssembly);
        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(79, link.firstJunctionIndexInSecond());
        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(DUP, link.svType());

        // DEL with inserted bases
        posJunction = new Junction(CHR_1, 500, POS_ORIENT);
        negJunction = new Junction(CHR_1, 1000, NEG_ORIENT);

        String insertedBases = "GGGGGGG";
        int insertedBaseLength = insertedBases.length();
        firstAssemblyBases = firstRefBases + insertedBases + firstExtensionBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(posJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        secondAssemblyBases = secondExtensionBases + insertedBases + secondRefBases;

        secondAssembly = new JunctionAssembly(negJunction, secondAssemblyBases.getBytes(), baseQuals, 80 + insertedBaseLength);

        link = assemblyLinker.tryAssemblyOverlap(secondAssembly, firstAssembly);
        assertNotNull(link);
        assertEquals(insertedBases, link.insertedBases());
        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(79, link.firstJunctionIndexInSecond());
        assertEquals(DEL, link.svType());
    }

    @Test
    public void testAssemblyInvertedSplits()
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
        String insertedBases = "GGGGGGG";
        int insertedBaseLength = insertedBases.length();
        firstAssemblyBases = firstRefBases + insertedBases + firstExtensionBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, 99);

        secondAssemblyBases = secondRefBases + reverseComplementBases(insertedBases) + secondExtensionBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, 99);

        link = assemblyLinker.tryAssemblyOverlap(firstAssembly, secondAssembly);
        assertNotNull(link);
        assertEquals(insertedBases, link.insertedBases());

        assertEquals(107, link.firstJunctionIndexInSecond());
        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(INV, link.svType());
    }
}
