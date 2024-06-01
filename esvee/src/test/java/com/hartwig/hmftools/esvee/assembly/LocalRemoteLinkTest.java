package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.formTestRefSequence;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;
import static com.hartwig.hmftools.esvee.assembly.phase.RemoteRegionAssembler.isExtensionCandidateAssembly;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.DISCORDANT;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.phase.LocalSequenceMatcher;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteRegionAssembler;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;

import org.junit.Test;

public class LocalRemoteLinkTest
{
    @Test
    public void testAssemblyLocalRefMatch()
    {
        // must be long enough to test the local ref genome sequence but not repetitive
        String localRefSequence = formTestRefSequence(600);

        MockRefGenome refGenome = new MockRefGenome(false);
        refGenome.RefGenomeMap.put(CHR_1, localRefSequence);

        LocalSequenceMatcher localSequenceMatcher = new LocalSequenceMatcher(refGenome, 200);

        // first a basic exact match junction
        Junction posJunction = new Junction(CHR_1, 300, FORWARD);

        String assemblyRefBases = localRefSequence.substring(201, 301);
        String assemblyExtensionBases = localRefSequence.substring(400, 450);
        String assemblyBases = assemblyRefBases + assemblyExtensionBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());

        JunctionAssembly assembly = new JunctionAssembly(posJunction, assemblyBases.getBytes(), baseQuals, assemblyRefBases.length() - 1);

        AssemblyLink assemblyLink = localSequenceMatcher.tryLocalAssemblyLink(assembly);
        assertNotNull(assemblyLink);
        assertEquals(DEL, assemblyLink.svType());
        assertEquals(400, assemblyLink.second().junction().Position);
        assertEquals(50, assemblyLink.second().refBaseLength());

        Junction negJunction = new Junction(CHR_1, 300, REVERSE);

        assemblyRefBases = localRefSequence.substring(300, 400);
        assemblyExtensionBases = localRefSequence.substring(400, 450);
        assemblyBases = assemblyExtensionBases + assemblyRefBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());

        assembly = new JunctionAssembly(negJunction, assemblyBases.getBytes(), baseQuals, assemblyExtensionBases.length());

        assemblyLink = localSequenceMatcher.tryLocalAssemblyLink(assembly);
        assertNotNull(assemblyLink);
        assertEquals(DUP, assemblyLink.svType());
        assertEquals(400, assemblyLink.second().junction().Position);
        assertEquals(50, assemblyLink.second().refBaseLength());
    }

    @Test
    public void testAssemblyRemoteReadMatch()
    {
        String remoteRefSequence = formTestRefSequence(400);

        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_2, remoteRefSequence);

        // first a basic exact match at the remote site
        String assemblyRefBases = REF_BASES_200.substring(1, 101);
        String assemblyExtensionBases = refGenome.getBaseString(CHR_2, 200, 300);
        String assemblyBases = assemblyRefBases + assemblyExtensionBases;

        JunctionAssembly assembly = createAssembly(CHR_1, 100, FORWARD, assemblyBases, 99);

        String juncReadBases = REF_BASES_200.substring(51, 101) + refGenome.getBaseString(CHR_2, 200, 250);

        Read juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 51, juncReadBases, "50M50S", CHR_2, 200, true);

        Read juncRead2 = cloneRead(juncRead, READ_ID_GENERATOR.nextId());

        assembly.addJunctionRead(juncRead);
        assembly.addJunctionRead(juncRead2);

        assertTrue(isExtensionCandidateAssembly(assembly));

        RemoteRegionAssembler remoteRegionAssembler = new RemoteRegionAssembler(refGenome, null);

        Read remoteRead = createRead(READ_ID_GENERATOR.nextId(), 200, refGenome.getBaseString(CHR_2, 200, 300), "100M");

        RemoteRegion remoteRegion = new RemoteRegion(
                new ChrBaseRegion(CHR_2, 200, 300), REVERSE, remoteRead.id(), DISCORDANT);

        remoteRegionAssembler.addMatchedReads(List.of(remoteRead), remoteRegion);

        byte[] remoteRegionBases = refGenome.getBases(CHR_2, 200, 300);

        AssemblyLink assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(
                assembly, 200, 300, remoteRegionBases);

        assertNotNull(assemblyLink);
        assertEquals(BND, assemblyLink.svType());
        assertEquals(200, assemblyLink.second().junction().Position);
        assertEquals(REVERSE, assemblyLink.second().junction().Orient);
        assertEquals(1, assemblyLink.second().supportCount());

        // now a match upstream of the remote read, inferring the missing ref bases
        assemblyExtensionBases = refGenome.getBaseString(CHR_2, 150, 250);
        assemblyBases = assemblyRefBases + assemblyExtensionBases;

        assembly = createAssembly(CHR_1, 100, FORWARD, assemblyBases, 99);

        assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(assembly, 200, 300, remoteRegionBases);

        assertNotNull(assemblyLink);
        assertEquals(BND, assemblyLink.svType());
        assertEquals(150, assemblyLink.second().junction().Position);
        assertEquals(REVERSE, assemblyLink.second().junction().Orient);
        assertEquals(1, assemblyLink.second().supportCount());

        // now test negative to negative orientation (a la 3-6-3)
        assemblyRefBases = REF_BASES_200.substring(100, 200);
        assemblyExtensionBases = refGenome.getBaseString(CHR_2, 150, 249);
        String assemblyExtensionBasesReversed = Nucleotides.reverseComplementBases(assemblyExtensionBases);
        assemblyBases = assemblyExtensionBasesReversed + assemblyRefBases;

        assembly = createAssembly(CHR_1, 100, REVERSE, assemblyBases, assemblyExtensionBasesReversed.length());

        assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(assembly, 200, 300, remoteRegionBases);

        assertNotNull(assemblyLink);
        assertEquals(BND, assemblyLink.svType());
        assertEquals(CHR_2, assemblyLink.first().junction().Chromosome);
        assertEquals(150, assemblyLink.first().junction().Position);
        assertEquals(REVERSE, assemblyLink.first().junction().Orient);
        assertEquals(1, assemblyLink.first().supportCount());

        assertEquals(100, assemblyLink.second().junction().Position);
        assertEquals(REVERSE, assemblyLink.second().junction().Orient);
    }
}
