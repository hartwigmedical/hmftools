package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.esvee.TestUtils.createConcordantRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.setMateCigar;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.JUNCTION_MATE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteGroupBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteReadExtractor;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.output.PhaseGroupBuildWriter;
import com.hartwig.hmftools.esvee.assembly.phase.LocalGroupBuilder;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.common.utils.TaskQueue;

import org.junit.Test;

public class PhaseGroupTest
{
    @Test
    public void testLocalPhaseGroupBuilding()
    {
        Junction posJunction1 = new Junction(CHR_1, 100, FORWARD);
        Junction negJunction1 = new Junction(CHR_1, 200, REVERSE);

        Junction posJunction2 = new Junction(CHR_1, 1500, FORWARD);
        Junction negJunction2 = new Junction(CHR_1, 2700, REVERSE);

        Junction posJunction3 = new Junction(CHR_1, 4000, FORWARD);
        Junction negJunction3 = new Junction(CHR_1, 4200, REVERSE);
        Junction negJunction4 = new Junction(CHR_1, 4500, REVERSE);
        Junction posJunction4 = new Junction(CHR_1, 4600, FORWARD);

        String assemblyBases = REF_BASES_200.substring(0, 100);

        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());

        JunctionAssembly posAssembly1 = new JunctionAssembly(posJunction1, assemblyBases.getBytes(), baseQuals, 50);

        Read read1 = createConcordantRead(READ_ID_GENERATOR.nextId(), 50, assemblyBases, "50M50S", 250);
        posAssembly1.addJunctionRead(read1);

        JunctionAssembly negAssembly1 = new JunctionAssembly(negJunction1, assemblyBases.getBytes(), baseQuals, 50);
        Read read2 = createConcordantRead(READ_ID_GENERATOR.nextId(), 200, assemblyBases, "50S50M", 10);
        setMateCigar(read2, "50M");
        negAssembly1.addJunctionRead(read2);

        JunctionAssembly posAssembly2 = new JunctionAssembly(posJunction2, assemblyBases.getBytes(), baseQuals, 50);
        JunctionAssembly negAssembly2 = new JunctionAssembly(negJunction2, assemblyBases.getBytes(), baseQuals, 50);

        JunctionAssembly posAssembly3 = new JunctionAssembly(posJunction3, assemblyBases.getBytes(), baseQuals, 50);

        Read read3 = createConcordantRead(READ_ID_GENERATOR.nextId(), 3950, assemblyBases, "50M50S", 4800);
        posAssembly3.addJunctionRead(read3);

        JunctionAssembly negAssembly3 = new JunctionAssembly(negJunction3, assemblyBases.getBytes(), baseQuals, 50);
        Read read4 = createConcordantRead(READ_ID_GENERATOR.nextId(), 4200, assemblyBases, "50S50M", 3500);
        setMateCigar(read4, "100M");
        negAssembly3.addJunctionRead(read4);

        JunctionAssembly posAssembly4 = new JunctionAssembly(posJunction4, assemblyBases.getBytes(), baseQuals, 50);

        Read read5 = createConcordantRead(READ_ID_GENERATOR.nextId(), 4550, assemblyBases, "50M50S", 3800);
        posAssembly4.addJunctionRead(read5);

        JunctionAssembly negAssembly4 = new JunctionAssembly(negJunction4, assemblyBases.getBytes(), baseQuals, 50);
        Read read6 = createConcordantRead(READ_ID_GENERATOR.nextId(), 4500, assemblyBases, "50S50M", 3800);
        setMateCigar(read6, "100M");
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

        LocalGroupBuilder builder = new LocalGroupBuilder(TEST_CONFIG, new TaskQueue(junctionGroups), writer);

        builder.run();

        Set<PhaseGroup> phaseGroups = builder.phaseGroups();
        assertEquals(2, phaseGroups.size());

        PhaseGroup phaseGroup = phaseGroups.stream().filter(x -> x.assemblies().contains(posAssembly1)).findFirst().orElse(null);
        assertNotNull(phaseGroup);
        assertTrue(phaseGroup.assemblies().contains(negAssembly1));

        phaseGroup = phaseGroups.stream().filter(x -> x.assemblies().contains(posAssembly3)).findFirst().orElse(null);
        assertNotNull(phaseGroup);
        assertEquals(4, phaseGroup.assemblyCount());
        assertTrue(phaseGroup.assemblies().contains(posAssembly3));
        assertTrue(phaseGroup.assemblies().contains(posAssembly4));
        assertTrue(phaseGroup.assemblies().contains(negAssembly3));
        assertTrue(phaseGroup.assemblies().contains(negAssembly4));
    }

    @Test
    public void testRemotePhaseGroupBuilding()
    {
        Queue<JunctionGroup> junctionGroups = new ConcurrentLinkedQueue<>();
        PhaseGroupBuildWriter writer = new PhaseGroupBuildWriter(TEST_CONFIG);

        Map<String,List<JunctionGroup>> junctionGroupMap = Maps.newHashMap();

        String assemblyBases = REF_BASES_200.substring(0, 50);
        JunctionAssembly assembly1 = createAssembly(CHR_1, 100, FORWARD, assemblyBases, 25);
        JunctionAssembly assembly2 = createAssembly(CHR_1, 500, FORWARD, assemblyBases, 25);
        JunctionAssembly assembly3 = createAssembly(CHR_1, 4100, FORWARD, assemblyBases, 25);
        JunctionAssembly assembly4 = createAssembly(CHR_1, 6100, FORWARD, assemblyBases, 25);
        JunctionAssembly assembly5 = createAssembly(CHR_1, 6500, FORWARD, assemblyBases, 25);

        JunctionGroup junctionGroup1 = new JunctionGroup(assembly1.junction());
        junctionGroup1.addJunction(assembly2.junction());
        junctionGroup1.addJunctionAssemblies(List.of(assembly1, assembly2));

        JunctionGroup junctionGroup2 = new JunctionGroup(assembly3.junction());
        junctionGroup2.addJunctionAssemblies(List.of(assembly3));

        JunctionGroup junctionGroup3 = new JunctionGroup(assembly4.junction());
        junctionGroup3.addJunction(assembly5.junction());
        junctionGroup3.addJunctionAssemblies(List.of(assembly4, assembly5));

        junctionGroupMap.put(CHR_1, List.of(junctionGroup1, junctionGroup2, junctionGroup3));

        JunctionAssembly assembly6 = createAssembly(CHR_2, 100, REVERSE, assemblyBases, 25);
        JunctionAssembly assembly7 = createAssembly(CHR_2, 500, REVERSE, assemblyBases, 25);
        JunctionAssembly assembly8 = createAssembly(CHR_2, 4100, REVERSE, assemblyBases, 25);
        JunctionAssembly assembly9 = createAssembly(CHR_2, 6100, REVERSE, assemblyBases, 25);
        JunctionAssembly assembly10 = createAssembly(CHR_2, 6500, REVERSE, assemblyBases, 25);

        JunctionGroup junctionGroup4 = new JunctionGroup(assembly6.junction());
        junctionGroup4.addJunction(assembly7.junction());
        junctionGroup4.addJunctionAssemblies(List.of(assembly6, assembly7));

        JunctionGroup junctionGroup5 = new JunctionGroup(assembly8.junction());
        junctionGroup5.addJunctionAssemblies(List.of(assembly1, assembly8));

        JunctionGroup junctionGroup6 = new JunctionGroup(assembly9.junction());
        junctionGroup6.addJunction(assembly10.junction());
        junctionGroup6.addJunctionAssemblies(List.of(assembly9, assembly10));

        junctionGroupMap.put(CHR_2, List.of(junctionGroup4, junctionGroup5, junctionGroup6));

        junctionGroupMap.values().stream().forEach(x -> junctionGroups.addAll(x));

        RemoteGroupBuilder remoteGroupBuilder = new RemoteGroupBuilder(TEST_CONFIG, new TaskQueue(junctionGroups), junctionGroupMap, writer);

        remoteGroupBuilder.run();

        Set<PhaseGroup> phaseGroups = remoteGroupBuilder.phaseGroups();

        assertEquals(0, phaseGroups.size());

        // now run again with links between the groups

        String readId1 = READ_ID_GENERATOR.nextId();
        String readId2 = READ_ID_GENERATOR.nextId();
        String readId3 = READ_ID_GENERATOR.nextId();

        // sufficient linking criteria
        RemoteRegion region1 = new RemoteRegion(
                new ChrBaseRegion(CHR_1, 6000, 6200), readId1, DISCORDANT);
        region1.addReadDetails(readId2, 6000, 6100, JUNCTION_MATE);

        RemoteRegion region1b = new RemoteRegion(
                new ChrBaseRegion(CHR_2, 6000, 6100), readId3, DISCORDANT);
        assembly1.addRemoteRegions(List.of(region1, region1b));

        assembly4.addJunctionRead(createRead(readId1, 6000, assemblyBases, "50M"));
        assembly4.addJunctionRead(createRead(readId2, 6000, assemblyBases, "50M"));

        // insufficient read count overlap
        assembly10.addJunctionRead(createRead(readId3, 6000, assemblyBases, "50M"));

        // link between 3 assemblies
        String readId4 = READ_ID_GENERATOR.nextId();
        String readId5 = READ_ID_GENERATOR.nextId();
        String readId6 = READ_ID_GENERATOR.nextId();
        String readId7 = READ_ID_GENERATOR.nextId();

        RemoteRegion region2 = new RemoteRegion(
                new ChrBaseRegion(CHR_2, 4000, 4200), readId4, DISCORDANT);
        region2.addReadDetails(readId5, 4000, 4100, JUNCTION_MATE);
        assembly3.addRemoteRegions(List.of(region2));

        assembly8.addJunctionRead(createRead(readId4, 6000, assemblyBases, "50M"));
        assembly8.addJunctionRead(createRead(readId5, 6000, assemblyBases, "50M"));

        RemoteRegion region3 = new RemoteRegion(
                new ChrBaseRegion(CHR_1, 450, 600), readId6, DISCORDANT);
        region3.addReadDetails(readId7, 450, 600, JUNCTION_MATE);
        assembly8.addRemoteRegions(List.of(region3));

        assembly2.addJunctionRead(createRead(readId6, 500, assemblyBases, "50M"));
        assembly2.addJunctionRead(createRead(readId7, 500, assemblyBases, "50M"));


        junctionGroupMap.values().stream().forEach(x -> junctionGroups.addAll(x));
        remoteGroupBuilder = new RemoteGroupBuilder(TEST_CONFIG, new TaskQueue(junctionGroups), junctionGroupMap, writer);
        remoteGroupBuilder.run();

        phaseGroups = remoteGroupBuilder.phaseGroups();

        assertEquals(2, phaseGroups.size());

        PhaseGroup phaseGroup = phaseGroups.stream().filter(x -> x.assemblies().contains(assembly1)).findFirst().orElse(null);
        assertNotNull(phaseGroup);
        assertTrue(phaseGroup.assemblies().contains(assembly4));

        phaseGroup = phaseGroups.stream().filter(x -> x.assemblies().contains(assembly3)).findFirst().orElse(null);
        assertNotNull(phaseGroup);
        assertEquals(3, phaseGroup.assemblyCount());
        assertTrue(phaseGroup.assemblies().contains(assembly2));
        assertTrue(phaseGroup.assemblies().contains(assembly8));
    }

    @Test
    public void testPhasetSetMerging()
    {
        // 2 pairs of links are merged
        String assemblyBases1a = REF_BASES_400.substring(121, 201) + REF_BASES_400.substring(250, 300);
        String assemblyBases1b = REF_BASES_400.substring(151, 201) + REF_BASES_400.substring(250, 350);

        JunctionAssembly assembly1a = createAssembly(CHR_1, 200, FORWARD, assemblyBases1a, 80);
        JunctionAssembly assembly1b = createAssembly(CHR_1, 250, REVERSE, assemblyBases1b, 50);

        Read juncRead1a = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 121, assemblyBases1a, "80M50S", CHR_1, 200, false);
        assembly1a.addJunctionRead(juncRead1a);

        Read juncRead1b = createRead(
                juncRead1a.id(), CHR_1, 250, assemblyBases1b, "50S100M", CHR_1, 200, false);
        assembly1b.addJunctionRead(juncRead1b);

        String assemblyBases2a = REF_BASES_400.substring(101, 201) + REF_BASES_400.substring(250, 300); // longer at start
        String assemblyBases2b = REF_BASES_400.substring(151, 201) + REF_BASES_400.substring(250, 330); // shorter at end
        JunctionAssembly assembly2a = createAssembly(CHR_2, 200, FORWARD, assemblyBases2a, 99);
        JunctionAssembly assembly2b = createAssembly(CHR_2, 250, REVERSE, assemblyBases2b, 50);

        Read juncRead2a = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 121, assemblyBases1a, "100M50S", CHR_2, 200, false);
        assembly2a.addJunctionRead(juncRead2a);

        Read juncRead2b = createRead(
                juncRead2a.id(), CHR_2, 250, assemblyBases1b, "50S80M", CHR_2, 200, false);
        assembly2b.addJunctionRead(juncRead2b);

        PhaseGroup phaseGroup = new PhaseGroup(assembly1a, assembly1b);
        phaseGroup.addAssembly(assembly2a);
        phaseGroup.addAssembly(assembly2b);

        RefGenomeInterface refGenome = new MockRefGenome();
        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(refGenome, new RemoteReadExtractor(null), phaseGroup);
        phaseSetBuilder.buildPhaseSets();

        assertEquals(2, phaseGroup.phaseSets().size());
        PhaseSet phaseSet1 = phaseGroup.phaseSets().get(0);
        PhaseSet phaseSet2 = phaseGroup.phaseSets().get(1);

        phaseGroup.finalisePhaseSetAlignments();

        assertFalse(phaseSet1.merged());
        assertTrue(phaseSet1.mergedPhaseSets().contains(phaseSet2));
        assertTrue(phaseSet2.merged());

        assertEquals(200, phaseSet1.assemblyAlignment().fullSequenceLength());

        // cannot merge if neither junction is overlapped - ie if just matching on ref bases
        assemblyBases1a = REF_BASES_400.substring(1, 51) + REF_BASES_400.substring(100, 150);
        assemblyBases1b = REF_BASES_400.substring(1, 51) + REF_BASES_400.substring(100, 250);

        assembly1a = createAssembly(CHR_1, 50, FORWARD, assemblyBases1a, 49);
        assembly1b = createAssembly(CHR_1, 100, REVERSE, assemblyBases1b, 50);

        juncRead1a = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, assemblyBases1a, "51M50S", CHR_1, 200, false);
        assembly1a.addJunctionRead(juncRead1a);

        juncRead1b = createRead(
                juncRead1a.id(), CHR_1, 100, assemblyBases1b, "50S151M", CHR_1, 50, false);
        assembly1b.addJunctionRead(juncRead1b);

        assemblyBases2a = REF_BASES_400.substring(120, 250) + REF_BASES_400.substring(300, 350);
        assemblyBases2b = REF_BASES_400.substring(200, 250) + REF_BASES_400.substring(300, 400);
        assembly2a = createAssembly(CHR_2, 250, FORWARD, assemblyBases2a, 149);
        assembly2b = createAssembly(CHR_2, 300, REVERSE, assemblyBases2b, 50);

        juncRead2a = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 120, assemblyBases1a, "131M50S", CHR_2, 200, false);
        assembly2a.addJunctionRead(juncRead2a);

        juncRead2b = createRead(
                juncRead2a.id(), CHR_2, 200, assemblyBases1b, "51S80M", CHR_2, 200, false);
        assembly2b.addJunctionRead(juncRead2b);

        phaseGroup = new PhaseGroup(assembly1a, assembly1b);
        phaseGroup.addAssembly(assembly2a);
        phaseGroup.addAssembly(assembly2b);

        phaseSetBuilder = new PhaseSetBuilder(refGenome, new RemoteReadExtractor(null), phaseGroup);
        phaseSetBuilder.buildPhaseSets();

        assertEquals(2, phaseGroup.phaseSets().size());
        phaseSet1 = phaseGroup.phaseSets().get(0);
        phaseSet2 = phaseGroup.phaseSets().get(1);

        phaseGroup.finalisePhaseSetAlignments();

        assertFalse(phaseSet1.merged());
        assertFalse(phaseSet2.merged());
    }
}
