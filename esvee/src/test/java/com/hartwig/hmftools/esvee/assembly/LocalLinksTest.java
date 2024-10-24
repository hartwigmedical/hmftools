package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_600;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.formTestRefSequence;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.AssemblyDeduper.dedupProximateAssemblies;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
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
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());

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
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());

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
        // local 32+ indel, junction matches with a nearby standard junction, is deduped and the other links
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_600);

        // the INDEL junction assemblies
        Junction posIndelJunction = new Junction(CHR_1, 100, FORWARD);
        posIndelJunction.markAsIndel();

        String refBases1 = refGenome.getBaseString(CHR_1, 1, 100);
        String duplicatedBases = refGenome.getBaseString(CHR_1, 101, 150);
        String extBases1 = duplicatedBases;
        String assemblyBases1 = refBases1 + extBases1;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases1.length());

        IndelCoords indelCoords = new IndelCoords(100, 101, 50);
        indelCoords.setInsertedBases(duplicatedBases);

        JunctionAssembly posIndelAssembly = new JunctionAssembly(posIndelJunction, assemblyBases1.getBytes(), baseQuals, refBases1.length() - 1);
        posIndelAssembly.setIndelCoords(indelCoords);

        Junction negIndelJunction = new Junction(CHR_1, 101, REVERSE);
        negIndelJunction.markAsIndel();

        String refBases2 = refGenome.getBaseString(CHR_1,101, 200);
        String extBases2 = refGenome.getBaseString(CHR_1,51, 100);
        String assemblyBases2 = extBases2 + refBases2;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases2.length());

        JunctionAssembly negIndelAssembly = new JunctionAssembly(negIndelJunction, assemblyBases2.getBytes(), baseQuals, extBases2.length());
        negIndelAssembly.setIndelCoords(indelCoords);

        // and the matching assembly for the
        Junction posJunction = new Junction(CHR_1, 150, FORWARD);

        String refBases3 = refGenome.getBaseString(CHR_1, 50, 150);
        String extBases3 = duplicatedBases;
        String assemblyBases3 = refBases3 + extBases3;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases1.length());

        JunctionAssembly posAssembly = new JunctionAssembly(posJunction, assemblyBases3.getBytes(), baseQuals, refBases3.length() - 1);

        Read juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 50, assemblyBases3,
                makeCigarString(refBases3, 0, extBases3.length()), CHR_1, 200, true);

        posAssembly.addJunctionRead(juncRead);

        Read juncMate = createRead(
                juncRead.id(), CHR_1, 101, assemblyBases2,
                makeCigarString(refBases2, extBases2.length(), 0), CHR_1, 200, false);
        juncMate.bamRecord().setReadNegativeStrandFlag(true);
        negIndelAssembly.addJunctionRead(juncMate);

        List<JunctionAssembly> junctionGroupAssemblies = Lists.newArrayList(posIndelAssembly, negIndelAssembly);
        List<JunctionAssembly> candidateAssemblies = Lists.newArrayList(posAssembly);
        dedupProximateAssemblies(junctionGroupAssemblies, candidateAssemblies);
        junctionGroupAssemblies.addAll(candidateAssemblies);

        assertEquals(2, junctionGroupAssemblies.size());
        assertTrue(junctionGroupAssemblies.contains(posAssembly));
        assertTrue(junctionGroupAssemblies.contains(negIndelAssembly));

        // now test linking logic
        PhaseGroup phaseGroup = new PhaseGroup(posAssembly, negIndelAssembly);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(
                refGenome, new RemoteRegionAssembler(refGenome, null), phaseGroup);

        phaseSetBuilder.buildPhaseSets();

        /* FIXME: change the ref sequence bases to look like a duplicated section
        assertEquals(1, phaseGroup.phaseSets().size());
        PhaseSet phaseSet = phaseGroup.phaseSets().get(0);

        assertEquals(2, phaseSet.assemblies().size());
        AssemblyLink link = phaseSet.assemblyLinks().get(0);

        assertNotNull(link);
        assertTrue(link.insertedBases().isEmpty());

        assertEquals(negIndelJunction, link.first());
        assertEquals(posAssembly, link.second());
        assertEquals(DUP, link.svType());
        */


    }
}
