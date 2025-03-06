package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_600;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.formTestRefSequence;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.AssemblyDeduper.dedupProximateAssemblies;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.tryAssemblyOverlap;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.FACING;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.SPLIT;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.phase.LocalSequenceMatcher;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetTask;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteRegionAssembler;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.common.FilterType;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import org.junit.Test;

public class LocalLinksTest
{
    @Test
    public void testAssemblyLocalRefMatch()
    {
        // must be long enough to test the local ref genome sequence but not repetitive
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_600);

        LocalSequenceMatcher localSequenceMatcher = new LocalSequenceMatcher(refGenome, 200);

        // first a basic exact match junction
        Junction posJunction = new Junction(CHR_1, 300, FORWARD);

        String assemblyRefBases = refGenome.getBaseString(CHR_1, 200, 300);
        String assemblyExtensionBases = refGenome.getBaseString(CHR_1, 400, 449);
        String assemblyBases = assemblyRefBases + assemblyExtensionBases;
        byte[] baseQuals = buildDefaultBaseQuals(assemblyBases.length());

        JunctionAssembly assembly = new JunctionAssembly(posJunction, assemblyBases.getBytes(), baseQuals, assemblyRefBases.length() - 1);

        AssemblyLink assemblyLink = localSequenceMatcher.tryLocalAssemblyLink(assembly);
        assertNotNull(assemblyLink);
        assertEquals(DEL, assemblyLink.svType());
        JunctionAssembly localAssembly = assemblyLink.second();

        assertEquals(400, localAssembly.junction().Position);
        assertTrue(localAssembly.junction().Orient.isReverse());
        assertEquals(50, localAssembly.extensionLength());
        String localSequence = refGenome.getBaseString(CHR_1, 251, 300) + refGenome.getBaseString(CHR_1, 400, 449);
        assertEquals(localSequence, localAssembly.formFullSequence());

        Junction negJunction = new Junction(CHR_1, 300, REVERSE);

        assemblyRefBases = refGenome.getBaseString(CHR_1,300, 400);
        assemblyExtensionBases = refGenome.getBaseString(CHR_1,401, 450);
        assemblyBases = assemblyExtensionBases + assemblyRefBases;
        baseQuals = buildDefaultBaseQuals(assemblyBases.length());

        assembly = new JunctionAssembly(negJunction, assemblyBases.getBytes(), baseQuals, assemblyExtensionBases.length());

        assemblyLink = localSequenceMatcher.tryLocalAssemblyLink(assembly);
        assertNotNull(assemblyLink);
        assertEquals(DUP, assemblyLink.svType());

        localAssembly = assemblyLink.first();
        assertEquals(450, localAssembly.junction().Position);
        assertTrue(localAssembly.junction().Orient.isForward());
        assertEquals(50, localAssembly.extensionLength());

        localSequence = refGenome.getBaseString(CHR_1,401, 450) + refGenome.getBaseString(CHR_1,300, 349);
        assertEquals(localSequence, localAssembly.formFullSequence());
    }

    @Test
    public void testLocalIndelVsJunction()
    {
        // 2 DELs with the internal long section having a duplication of 35 bases, represented as both an INDEL and a soft-clip assembly
        // first DEL:   1:100-200 - assemblies 1 & 2
        // INDEL:       1:250:251 - assemblies 3 & 4, with insert of 35 bases matching the sequence 251-285
        // dedup asm    1:285:1 - assembly 5
        // second DEL:  1:350-450 - assemblies 6 & 7
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_600);

        Junction junction1 = new Junction(CHR_1, 100, FORWARD);

        String refBases1 = refGenome.getBaseString(CHR_1, 1, 100);
        String extBases1 = refGenome.getBaseString(CHR_1, 200, 249);
        String assemblyBases1 = refBases1 + extBases1;

        JunctionAssembly assembly1 = new JunctionAssembly(
                junction1, assemblyBases1.getBytes(), buildDefaultBaseQuals(assemblyBases1.length()), refBases1.length() - 1);

        Junction junction2 = new Junction(CHR_1, 200, REVERSE);

        String refBases2 = refGenome.getBaseString(CHR_1, 200, 250);
        String extBases2 = refGenome.getBaseString(CHR_1, 51, 100);
        String assemblyBases2 = extBases2 + refBases2;

        JunctionAssembly assembly2 = new JunctionAssembly(
                junction2, assemblyBases2.getBytes(), buildDefaultBaseQuals(assemblyBases2.length()), extBases2.length());

        // read pair linking assemblies 1 & 2
        Read juncRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, assemblyBases1,
                makeCigarString(refBases1, 0, extBases1.length()), CHR_1, 200, true);

        assembly1.addJunctionRead(juncRead1);

        Read juncRead2 = createRead(
                juncRead1.id(), CHR_1, 200, assemblyBases2,
                makeCigarString(refBases2, extBases2.length(), 0), CHR_1, 1, false);

        assembly2.addJunctionRead(juncRead2);

        Junction junction3 = new Junction(CHR_1, 250, FORWARD);
        junction3.markAsIndel();

        String duplicatedBases = refGenome.getBaseString(CHR_1, 251, 285);
        IndelCoords indelCoords = new IndelCoords(250, 251, duplicatedBases.length());
        indelCoords.setInsertedBases(duplicatedBases);

        String refBases3 = refGenome.getBaseString(CHR_1, 200, 250);
        String extBases3 = duplicatedBases;
        String assemblyBases3 = refBases3 + extBases3;

        JunctionAssembly assembly3 = new JunctionAssembly(
                junction3, assemblyBases3.getBytes(), buildDefaultBaseQuals(assemblyBases3.length()), refBases3.length() - 1);
        assembly3.setIndelCoords(indelCoords);

        Junction junction4 = new Junction(CHR_1, 251, REVERSE);
        junction4.markAsIndel();

        String refBases4 = refGenome.getBaseString(CHR_1, 251, 350);
        String extBases4 = duplicatedBases;
        String assemblyBases4 = extBases4 + refBases4;

        JunctionAssembly assembly4 = new JunctionAssembly(
                junction4, assemblyBases4.getBytes(), buildDefaultBaseQuals(assemblyBases4.length()), extBases4.length());
        assembly4.setIndelCoords(indelCoords);

        // indel read and read linking to assembly 2
        Read juncRead2b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 200, assemblyBases2,
                makeCigarString(refBases2, extBases2.length(), 0), CHR_1, 221, true);

        assembly2.addJunctionRead(juncRead2b);

        String read3Bases = refGenome.getBaseString(CHR_1, 221, 250) + duplicatedBases
                + refGenome.getBaseString(CHR_1, 251, 280);

        Read juncRead3 = createRead(
                juncRead2b.id(), CHR_1, 221, read3Bases,
                "30M35I30M", CHR_1, 200, false);

        assembly3.addJunctionRead(juncRead3);
        assembly4.addJunctionRead(juncRead3);

        Junction junction5 = new Junction(CHR_1, 285, FORWARD);

        String refBases5 = refGenome.getBaseString(CHR_1, 200, 285);
        String extBases5 = duplicatedBases;
        String assemblyBases5 = refBases5 + extBases5;

        JunctionAssembly assembly5 = new JunctionAssembly(
                junction5, assemblyBases5.getBytes(), buildDefaultBaseQuals(assemblyBases5.length()), refBases5.length() - 1);

        Junction junction6 = new Junction(CHR_1, 350, FORWARD);

        String refBases6 = refGenome.getBaseString(CHR_1, 251, 350);
        String extBases6 = refGenome.getBaseString(CHR_1, 450, 500);
        String assemblyBases6 = refBases6 + extBases6;

        JunctionAssembly assembly6 = new JunctionAssembly(
                junction6, assemblyBases6.getBytes(), buildDefaultBaseQuals(assemblyBases6.length()), refBases6.length() - 1);

        String read4Bases = refGenome.getBaseString(CHR_1, 221, 250) + duplicatedBases
                + refGenome.getBaseString(CHR_1, 251, 280);

        Read juncRead4 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 221, read4Bases,
                "30M35I30M", CHR_1, 200, true);

        assembly3.addJunctionRead(juncRead4);
        assembly4.addJunctionRead(juncRead4);

        Read juncRead6 = createRead(
                juncRead4.id(), CHR_1, 251, assemblyBases6,
                makeCigarString(refBases6, 0, extBases6.length()), CHR_1, 221, false);

        assembly6.addJunctionRead(juncRead6);

        Junction junction7 = new Junction(CHR_1, 450, REVERSE);

        String refBases7 = refGenome.getBaseString(CHR_1, 450, 550);
        String extBases7 = refGenome.getBaseString(CHR_1, 300, 350);
        String assemblyBases7 = extBases7 + refBases7;

        JunctionAssembly assembly7 = new JunctionAssembly(
                junction7, assemblyBases7.getBytes(), buildDefaultBaseQuals(assemblyBases7.length()), extBases7.length());

        // read pair linking assemblies 6 & 7
        Read juncRead6b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 251, assemblyBases6,
                makeCigarString(refBases6, 0, extBases6.length()), CHR_1, 450, true);

        assembly1.addJunctionRead(juncRead1);

        Read juncRead7 = createRead(
                juncRead6b.id(), CHR_1, 450, assemblyBases2,
                makeCigarString(refBases7, extBases7.length(), 0), CHR_1, 251, false);

        assembly7.addJunctionRead(juncRead7);

        List<JunctionAssembly> junctionGroupAssemblies = Lists.newArrayList(assembly3, assembly4);
        List<JunctionAssembly> candidateAssemblies = Lists.newArrayList(assembly5);
        dedupProximateAssemblies(junctionGroupAssemblies, candidateAssemblies);
        junctionGroupAssemblies.addAll(candidateAssemblies);

        assertEquals(2, junctionGroupAssemblies.size());
        assertTrue(junctionGroupAssemblies.contains(assembly3));
        assertTrue(junctionGroupAssemblies.contains(assembly4));

        // now test linking logic
        PhaseGroup phaseGroup = new PhaseGroup(assembly1, assembly2);
        phaseGroup.addAssembly(assembly3);
        phaseGroup.addAssembly(assembly4);
        phaseGroup.addAssembly(assembly6);
        phaseGroup.addAssembly(assembly7);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(
                refGenome, new RemoteRegionAssembler(refGenome, null), phaseGroup);

        phaseSetBuilder.buildPhaseSets();

        assertEquals(1, phaseGroup.phaseSets().size());
        PhaseSet phaseSet = phaseGroup.phaseSets().get(0);

        assertEquals(6, phaseSet.assemblies().size());
        AssemblyLink link = phaseSet.assemblyLinks().get(0);

        assertEquals(SPLIT, link.type());
        assertEquals(assembly1, link.first());
        assertEquals(assembly2, link.second());

        link = phaseSet.assemblyLinks().get(1);
        assertEquals(FACING, link.type());
        assertEquals(assembly2, link.first());
        assertEquals(assembly3, link.second());

        link = phaseSet.assemblyLinks().get(2);
        assertEquals(INDEL, link.type());
        assertEquals(assembly3, link.first());
        assertEquals(assembly4, link.second());
        assertEquals(duplicatedBases, link.insertedBases());

        link = phaseSet.assemblyLinks().get(3);
        assertEquals(FACING, link.type());
        assertEquals(assembly4, link.first());
        assertEquals(assembly6, link.second());

        link = phaseSet.assemblyLinks().get(4);
        assertEquals(SPLIT, link.type());
        assertEquals(assembly6, link.first());
        assertEquals(assembly7, link.second());
    }

    @Test
    public void testLocalIndelRequiredOverlaps()
    {
        // test 1: there is sufficient overlap in ref and extension bases to form a link
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_600);

        // modelled as a DEL 1:100:1 to 1:300:-1
        Junction junction1 = new Junction(CHR_1, 100, FORWARD);

        String refBases1 = refGenome.getBaseString(CHR_1, 1, 100);
        String extBases1 = refGenome.getBaseString(CHR_1, 300, 349);
        String assemblyBases1 = refBases1 + extBases1;

        JunctionAssembly assembly1 = new JunctionAssembly(
                junction1, assemblyBases1.getBytes(), buildDefaultBaseQuals(assemblyBases1.length()), refBases1.length() - 1);

        Junction junction2 = new Junction(CHR_1, 300, REVERSE);

        String refBases2 = refGenome.getBaseString(CHR_1, 300, 399);
        String extBases2 = refGenome.getBaseString(CHR_1, 51, 100);
        String assemblyBases2 = extBases2 + refBases2;

        JunctionAssembly assembly2 = new JunctionAssembly(
                junction2, assemblyBases2.getBytes(), buildDefaultBaseQuals(assemblyBases2.length()), extBases2.length());

        AssemblyLink assemblyLink = tryAssemblyOverlap(assembly1, assembly2);
        assertNotNull(assemblyLink);

        // test 2: with a single mismatch, trigger the seed-search logic
        String extBases2Mod = extBases2.substring(0, 20) + getNextBase(extBases2.substring(20, 21)) + extBases2.substring(21);
        assemblyBases2 = extBases2Mod + refBases2;

        assembly2 = new JunctionAssembly(
                junction2, assemblyBases2.getBytes(), buildDefaultBaseQuals(assemblyBases2.length()), extBases2.length());

        assemblyLink = tryAssemblyOverlap(assembly1, assembly2);
        assertNotNull(assemblyLink);

        // test 3: sufficient overlap in ref and extension bases but too many mis-matches
        extBases2Mod = extBases2.substring(0, 20) + "TTTTT" + extBases2.substring(25);
        assemblyBases2 = extBases2Mod + refBases2;

        assembly2 = new JunctionAssembly(
                junction2, assemblyBases2.getBytes(), buildDefaultBaseQuals(assemblyBases2.length()), extBases2.length());

        assemblyLink = tryAssemblyOverlap(assembly1, assembly2);
        assertNull(assemblyLink);

        // test 4: insufficient overlap in ref and extension bases to form a link
        extBases2Mod = extBases2.substring(35);
        assemblyBases2 = extBases2Mod + refBases2;

        assembly2 = new JunctionAssembly(
                junction2, assemblyBases2.getBytes(), buildDefaultBaseQuals(assemblyBases2.length()), extBases2Mod.length());

        String refBases1Mod = extBases2.substring(35, 50) + refBases2.substring(0, 20);
        String extBases1Mod = refBases2.substring(20, 34) + getNextBase(refBases2.substring(34, 35));
        assemblyBases1 = refBases1Mod + extBases1Mod;

        assembly1 = new JunctionAssembly(
                junction1, assemblyBases1.getBytes(), buildDefaultBaseQuals(assemblyBases1.length()), refBases1Mod.length() - 1);

        assemblyLink = tryAssemblyOverlap(assembly1, assembly2);
        assertNull(assemblyLink);
    }
}
