package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
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
import com.hartwig.hmftools.esvee.assembly.phase.RemoteRegionAssembler;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;

import org.junit.Test;

public class RemoteLinksTest
{
    @Test
    public void testAssemblyRemoteReadMatch()
    {
        MockRefGenome refGenome = new MockRefGenome(false);
        refGenome.RefGenomeMap.put(CHR_2, REF_BASES_400);

        int remoteRegionStart = 200;
        int remoteRegionEnd = 300;
        String remoteRefBases = refGenome.getBaseString(CHR_2, remoteRegionStart, remoteRegionEnd);

        // test 1: extension bases match fully within the remote site - test each assembly and remote orientation combination in turn
        String assemblyRefBases = REF_BASES_200.substring(1, 101);
        String assemblyExtBases = refGenome.getBaseString(CHR_2, 220, 280);
        String assemblyBases = assemblyRefBases + assemblyExtBases;

        int juncPosition = 100;
        JunctionAssembly assembly = createAssembly(CHR_1, juncPosition, FORWARD, assemblyBases, assemblyRefBases.length() - 1);

        String juncReadBases = REF_BASES_200.substring(51, juncPosition + 1) + assemblyExtBases;

        // add reads to test out the criteria for attempting remote ref matching, does not factor into matching
        Read juncRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 51, juncReadBases,
                makeCigarString(juncReadBases, 0, assemblyExtBases.length()), CHR_2, 220, true);

        Read juncRead2 = cloneRead(juncRead, READ_ID_GENERATOR.nextId());

        assembly.addJunctionRead(juncRead);
        assembly.addJunctionRead(juncRead2);

        assertTrue(isExtensionCandidateAssembly(assembly));

        RemoteRegionAssembler remoteRegionAssembler = new RemoteRegionAssembler(refGenome, null);

        Read remoteRead = createRead(READ_ID_GENERATOR.nextId(), remoteRegionStart, remoteRefBases, makeCigarString(remoteRefBases, 0, 0));

        RemoteRegion remoteRegion = new RemoteRegion(
                new ChrBaseRegion(CHR_2, remoteRegionStart, remoteRegionEnd), remoteRead.id(), DISCORDANT);

        remoteRegionAssembler.addMatchedReads(List.of(remoteRead), remoteRegion);

        AssemblyLink assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(
                assembly, remoteRegion.start(), remoteRegion.end(), remoteRefBases.getBytes());

        assertNotNull(assemblyLink);
        assertEquals(BND, assemblyLink.svType());
        JunctionAssembly remoteAssembly = assemblyLink.otherAssembly(assembly);

        assertEquals(220, remoteAssembly.junction().Position);
        assertEquals(REVERSE, remoteAssembly.junction().Orient);
        assertEquals(1, remoteAssembly.supportCount());

        // test 2: assembly has reversed bases, requiring a switched orientation match in the remote region
        assemblyBases = assemblyRefBases + Nucleotides.reverseComplementBases(assemblyExtBases);

        assembly = createAssembly(CHR_1, juncPosition, FORWARD, assemblyBases, assemblyRefBases.length() - 1);

        assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(
                assembly, remoteRegion.start(), remoteRegion.end(), remoteRefBases.getBytes());

        assertNotNull(assemblyLink);
        remoteAssembly = assemblyLink.otherAssembly(assembly);

        assertEquals(280, remoteAssembly.junction().Position);
        assertEquals(FORWARD, remoteAssembly.junction().Orient);
        assertEquals(1, remoteAssembly.supportCount());

        // test 3: assembly is reversed, remote bases are not
        assemblyBases = assemblyExtBases + assemblyRefBases;

        assembly = createAssembly(CHR_1, juncPosition, REVERSE, assemblyBases, assemblyExtBases.length());

        assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(
                assembly, remoteRegion.start(), remoteRegion.end(), remoteRefBases.getBytes());

        assertNotNull(assemblyLink);
        remoteAssembly = assemblyLink.otherAssembly(assembly);

        assertEquals(280, remoteAssembly.junction().Position);
        assertEquals(FORWARD, remoteAssembly.junction().Orient);
        assertEquals(1, remoteAssembly.supportCount());

        // test 4: assembly is reversed as are the remote bases
        assemblyBases = Nucleotides.reverseComplementBases(assemblyExtBases) + assemblyRefBases;

        assembly = createAssembly(CHR_1, juncPosition, REVERSE, assemblyBases, assemblyExtBases.length());

        assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(
                assembly, remoteRegion.start(), remoteRegion.end(), remoteRefBases.getBytes());

        assertNotNull(assemblyLink);
        remoteAssembly = assemblyLink.otherAssembly(assembly);

        assertEquals(220, remoteAssembly.junction().Position);
        assertEquals(REVERSE, remoteAssembly.junction().Orient);
        assertEquals(1, remoteAssembly.supportCount());

        // test 5: extension bases match downstream of the remote region, but by extending it show an exact match
        assemblyExtBases = refGenome.getBaseString(CHR_2, 180, 280);
        assemblyBases = assemblyRefBases + assemblyExtBases;

        assembly = createAssembly(CHR_1, juncPosition, FORWARD, assemblyBases, assemblyRefBases.length() - 1);

        assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(
                assembly, remoteRegion.start(), remoteRegion.end(), remoteRefBases.getBytes());

        assertNotNull(assemblyLink);
        remoteAssembly = assemblyLink.otherAssembly(assembly);

        assertEquals(180, remoteAssembly.junction().Position);
        assertEquals(REVERSE, remoteAssembly.junction().Orient);
        assertEquals(1, remoteAssembly.supportCount());

        // test 6: extension bases match upstream of the remote region, but by extending it show an exact match
        assemblyExtBases = refGenome.getBaseString(CHR_2, 230, 350);
        assemblyBases = assemblyRefBases + assemblyExtBases;

        assembly = createAssembly(CHR_1, juncPosition, FORWARD, assemblyBases, assemblyRefBases.length() - 1);

        assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(
                assembly, remoteRegion.start(), remoteRegion.end(), remoteRefBases.getBytes());

        assertNotNull(assemblyLink);
        remoteAssembly = assemblyLink.otherAssembly(assembly);

        assertEquals(230, remoteAssembly.junction().Position);
        assertEquals(REVERSE, remoteAssembly.junction().Orient);
        assertEquals(1, remoteAssembly.supportCount());

        // test again with reversed orientations
        assemblyBases = assemblyExtBases + assemblyRefBases;

        assembly = createAssembly(CHR_1, juncPosition, REVERSE, assemblyBases, assemblyExtBases.length());

        assemblyLink = remoteRegionAssembler.tryAssemblyRemoteRefOverlap(
                assembly, remoteRegion.start(), remoteRegion.end(), remoteRefBases.getBytes());

        assertNotNull(assemblyLink);
        remoteAssembly = assemblyLink.otherAssembly(assembly);

        assertEquals(350, remoteAssembly.junction().Position);
        assertEquals(FORWARD, remoteAssembly.junction().Orient);
        assertEquals(1, remoteAssembly.supportCount());
    }
}
