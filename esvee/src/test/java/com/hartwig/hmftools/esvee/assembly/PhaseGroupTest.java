package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.esvee.TestUtils.createConcordantRead;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.output.PhaseGroupBuildWriter;
import com.hartwig.hmftools.esvee.assembly.phase.LocalGroupBuilder;
import com.hartwig.hmftools.esvee.assembly.read.Read;

import org.junit.Test;

public class PhaseGroupTest
{
    @Test
    public void testLocalPhaseGroupBuilding()
    {
        Junction posJunction1 = new Junction(CHR_1, 100, POS_ORIENT);
        Junction negJunction1 = new Junction(CHR_1, 200, NEG_ORIENT);

        Junction posJunction2 = new Junction(CHR_1, 1500, POS_ORIENT);
        Junction negJunction2 = new Junction(CHR_1, 2700, NEG_ORIENT);

        Junction posJunction3 = new Junction(CHR_1, 4000, POS_ORIENT);
        Junction negJunction3 = new Junction(CHR_1, 4200, NEG_ORIENT);
        Junction negJunction4 = new Junction(CHR_1, 4500, NEG_ORIENT);
        Junction posJunction4 = new Junction(CHR_1, 4600, POS_ORIENT);

        String assemblyBases = REF_BASES_200.substring(0, 100);

        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());

        JunctionAssembly posAssembly1 = new JunctionAssembly(posJunction1, assemblyBases.getBytes(), baseQuals, 50);

        Read read1 = createConcordantRead(READ_ID_GENERATOR.nextId(), 50, assemblyBases, "50M50S", 250);
        posAssembly1.addJunctionRead(read1);

        JunctionAssembly negAssembly1 = new JunctionAssembly(negJunction1, assemblyBases.getBytes(), baseQuals, 50);
        Read read2 = createConcordantRead(READ_ID_GENERATOR.nextId(), 200, assemblyBases, "50S50M", 10);
        negAssembly1.addJunctionRead(read2);

        JunctionAssembly posAssembly2 = new JunctionAssembly(posJunction2, assemblyBases.getBytes(), baseQuals, 50);
        JunctionAssembly negAssembly2 = new JunctionAssembly(negJunction2, assemblyBases.getBytes(), baseQuals, 50);

        JunctionAssembly posAssembly3 = new JunctionAssembly(posJunction3, assemblyBases.getBytes(), baseQuals, 50);

        Read read3 = createConcordantRead(READ_ID_GENERATOR.nextId(), 3950, assemblyBases, "50M50S", 4800);
        posAssembly3.addJunctionRead(read3);

        JunctionAssembly negAssembly3 = new JunctionAssembly(negJunction3, assemblyBases.getBytes(), baseQuals, 50);
        Read read4 = createConcordantRead(READ_ID_GENERATOR.nextId(), 4200, assemblyBases, "50S50M", 3500);
        negAssembly3.addJunctionRead(read4);

        JunctionAssembly posAssembly4 = new JunctionAssembly(posJunction4, assemblyBases.getBytes(), baseQuals, 50);

        Read read5 = createConcordantRead(READ_ID_GENERATOR.nextId(), 4550, assemblyBases, "50M50S", 4200);
        posAssembly4.addJunctionRead(read5);

        JunctionAssembly negAssembly4 = new JunctionAssembly(negJunction4, assemblyBases.getBytes(), baseQuals, 50);
        Read read6 = createConcordantRead(READ_ID_GENERATOR.nextId(), 4500, assemblyBases, "50S50M", 4200);
        negAssembly4.addJunctionRead(read6);

        Queue<JunctionGroup> junctionGroups = new ConcurrentLinkedQueue<>();
        PhaseGroupBuildWriter writer = new PhaseGroupBuildWriter(TEST_CONFIG);

        JunctionGroup junctionGroup = new JunctionGroup(posJunction1);
        junctionGroup.addJunction(posJunction2);
        junctionGroup.addJunction(posJunction3);
        junctionGroup.addJunction(posJunction4);
        junctionGroup.addJunction(negJunction1);
        junctionGroup.addJunction(negJunction2);
        junctionGroup.addJunction(negJunction3);
        junctionGroup.addJunction(negJunction4);

        junctionGroup.addJunctionAssemblies(List.of(
                posAssembly1, posAssembly2, posAssembly3, posAssembly4,
                negAssembly1, negAssembly2, negAssembly3, negAssembly4));

        junctionGroups.add(junctionGroup);

        LocalGroupBuilder builder = new LocalGroupBuilder(TEST_CONFIG, junctionGroups, writer);

        builder.run();

        Set<PhaseGroup> phaseGroups = builder.phaseGroups();
        assertEquals(2, phaseGroups.size());

        PhaseGroup phaseGroup = phaseGroups.stream().filter(x -> x.assemblies().contains(posAssembly1)).findFirst().orElse(null);
        assertNotNull(phaseGroup);
        assertTrue(phaseGroup.assemblies().contains(negAssembly1));

        phaseGroup = phaseGroups.stream().filter(x -> x.assemblies().contains(posAssembly3)).findFirst().orElse(null);
        assertNotNull(phaseGroup);
        assertEquals(3, phaseGroup.assemblyCount());
        assertTrue(phaseGroup.assemblies().contains(negAssembly3));
        assertTrue(phaseGroup.assemblies().contains(negAssembly4));
    }
}
