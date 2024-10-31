package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_600;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.getSupportTypeCount;
import static com.hartwig.hmftools.esvee.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.hasAssemblyLink;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.isValidSupportCoordsVsJunction;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.EXTENSION;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION_MATE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteRegionAssembler;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.RemoteReadType;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import org.junit.Test;

public class SpecificAlignmentsTest
{
    @Test
    public void testUnmappedSglPair()
    {
        // two facing SGLs which both have unmapped reads which form consistent extension sequences

        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_400);

        String refBases1 = refGenome.getBaseString(CHR_1, 200, 300);
        String extBases1 = REF_BASES_600.substring(420, 500);
        String assemblyBases1 = extBases1 + refBases1;

        JunctionAssembly assembly1 = createAssembly(CHR_1, 200, REVERSE, assemblyBases1, extBases1.length());

        String refBases2 = refGenome.getBaseString(CHR_1, 200, 300);
        String extBases2 = REF_BASES_600.substring(100, 180);
        String assemblyBases2 = refBases2 + extBases2;

        JunctionAssembly assembly2 = createAssembly(CHR_1, 300, FORWARD, assemblyBases2, refBases2.length() - 1);

        List<Read> candidates1 = Lists.newArrayList();
        List<Read> candidates2 = Lists.newArrayList();

        // a pair of mate junction reads
        Read juncRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 200, assemblyBases1.substring(0, 100), "50S50M",
                CHR_1, 251, true);
        juncRead1.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, "50M50S");

        Read juncRead2 = createRead(
                juncRead1.id(), CHR_1, 251, assemblyBases2.substring(51), "50M50S",
                CHR_1, 200, false);
        juncRead2.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, "50S50M");
        juncRead2.bamRecord().setSecondOfPairFlag(true);
        juncRead2.bamRecord().setFirstOfPairFlag(false);

        assembly1.addJunctionRead(juncRead1);
        assembly2.addJunctionRead(juncRead2);
        juncRead1.setMateRead(juncRead2);
        juncRead2.setMateRead(juncRead1);
        setSecondInPair(juncRead2.bamRecord());
        candidates1.add(juncRead2);
        candidates2.add(juncRead1);

        // a shared junction read with mate unmapped on the first assembly side
        String readBases = extBases1.substring(40) + refBases1 + extBases1.substring(0, 10);
        Read sharedRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 200, readBases, "10S101M10S",
                NO_CHROMOSOME_NAME, NO_POSITION, true);
        sharedRead1.bamRecord().setMateUnmappedFlag(true);

        assembly1.addJunctionRead(sharedRead1);
        assembly2.addJunctionRead(sharedRead1);

        // more junction reads whose mates extend out the unmapped sequence
        Read juncRead1b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 200, assemblyBases1.substring(0, 100), "50S50M",
                NO_CHROMOSOME_NAME, NO_POSITION, true);
        juncRead1b.bamRecord().setMateUnmappedFlag(true);
        juncRead1b.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, NO_CIGAR);

        readBases = REF_BASES_600.substring(400, 480); // extends the extension seq 20 bases further back
        Read unmappedRead1 = createRead(
                juncRead1b.id(), CHR_1, 200, readBases, NO_CIGAR, CHR_1, 200, false);
        unmappedRead1.bamRecord().setReadUnmappedFlag(true);
        unmappedRead1.bamRecord().setReadNegativeStrandFlag(true);

        assembly1.addJunctionRead(juncRead1b);
        candidates1.add(unmappedRead1);
        juncRead1b.setMateRead(unmappedRead1);
        unmappedRead1.setMateRead(juncRead1b);
        setSecondInPair(unmappedRead1.bamRecord());

        Read juncRead2b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 251, assemblyBases2.substring(51), "50M50S",
                NO_CHROMOSOME_NAME, NO_POSITION, true);
        juncRead2b.bamRecord().setMateUnmappedFlag(true);
        juncRead2b.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, NO_CIGAR);

        readBases = REF_BASES_600.substring(120, 200); // was 100-180, so extends 20 bases further out
        readBases = Nucleotides.reverseComplementBases(readBases);
        Read unmappedRead2 = createRead(
                juncRead2b.id(), CHR_1, 251, readBases, NO_CIGAR, CHR_1, 251, false);
        unmappedRead2.bamRecord().setReadUnmappedFlag(true);

        assembly2.addJunctionRead(juncRead2b);
        candidates2.add(unmappedRead2);
        juncRead2b.setMateRead(unmappedRead2);
        unmappedRead2.setMateRead(juncRead2b);
        setSecondInPair(unmappedRead2.bamRecord());


        RefBaseExtender refBaseExtender = new RefBaseExtender();
        refBaseExtender.findAssemblyCandidateExtensions(assembly1, candidates1);
        refBaseExtender.findAssemblyCandidateExtensions(assembly2, candidates2);

        PhaseGroup phaseGroup = new PhaseGroup(assembly2, assembly1);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(refGenome, new RemoteRegionAssembler(refGenome, null), phaseGroup);
        phaseSetBuilder.buildPhaseSets();

        // first check unmapped extensions
        assertEquals(5, assembly1.supportCount());
        assertEquals(4, assembly2.supportCount());

        assertEquals(100, assembly1.extensionLength());
        assertEquals(REF_BASES_600.substring(400, 500), assembly1.formJunctionSequence());
        assertEquals(100, assembly2.extensionLength());
        assertEquals(REF_BASES_600.substring(100, 200), assembly2.formJunctionSequence());

        SupportRead supportRead1 = assembly1.support().stream()
                .filter(x -> x.type() == EXTENSION && x.id().equals(unmappedRead1.id())).findFirst().orElse(null);
        assertNotNull(supportRead1);
        assertEquals(100, supportRead1.junctionReadStartDistance());

        SupportRead supportRead2 = assembly2.support().stream()
                .filter(x -> x.type() == EXTENSION && x.id().equals(unmappedRead2.id())).findFirst().orElse(null);
        assertNotNull(supportRead2);
        assertEquals(-20, supportRead2.junctionReadStartDistance());

        assertEquals(1, phaseGroup.phaseSets().size());
        PhaseSet phaseSet = phaseGroup.phaseSets().get(0);

        assertEquals(1, phaseSet.assemblyLinks().size());
        assertEquals(2, phaseSet.assemblies().size());

        hasAssemblyLink(phaseSet.assemblyLinks(), assembly1, assembly2, LinkType.FACING);

        AssemblyAlignment assemblyAlignment = new AssemblyAlignment(phaseSet);

        assertEquals(301, assemblyAlignment.fullSequenceLength());

        String fullSequence = REF_BASES_600.substring(400, 500) + refBases1 + REF_BASES_600.substring(100, 200);
        assertEquals(fullSequence, assemblyAlignment.fullSequence());

        // confirm full seq read indices
        SupportRead supportRead = assembly1.support().stream()
                .filter(x -> x.type() == JUNCTION && x.id().equals(juncRead1.id())).findFirst().orElse(null);
        assertNotNull(supportRead);
        assertEquals(50, supportRead.fullAssemblyIndexStart());

        supportRead = assembly1.support().stream()
                .filter(x -> x.type() == JUNCTION && x.id().equals(sharedRead1.id())).findFirst().orElse(null);
        assertNotNull(supportRead);
        assertEquals(90, supportRead.fullAssemblyIndexStart());

        supportRead = assembly1.support().stream()
                .filter(x -> x.type() == JUNCTION_MATE && x.id().equals(juncRead1.id())).findFirst().orElse(null);
        assertNotNull(supportRead);
        assertEquals(151, supportRead.fullAssemblyIndexStart());

        supportRead = assembly2.support().stream()
                .filter(x -> x.type() == JUNCTION && x.id().equals(juncRead1.id())).findFirst().orElse(null);
        assertNotNull(supportRead);
        assertEquals(151, supportRead.fullAssemblyIndexStart());

        supportRead = assembly2.support().stream()
                .filter(x -> x.type() == JUNCTION && x.id().equals(sharedRead1.id())).findFirst().orElse(null);
        assertNotNull(supportRead);
        assertEquals(90, supportRead.fullAssemblyIndexStart());

        supportRead = assembly2.support().stream()
                .filter(x -> x.type() == EXTENSION && x.id().equals(unmappedRead2.id())).findFirst().orElse(null);
        assertNotNull(supportRead);
        assertEquals(221, supportRead.fullAssemblyIndexStart());
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

        String refBases1 = refGenome.getBaseString(CHR_1, 1, 100);
        String extBases1 = refGenome.getBaseString(CHR_2, 50, 150);
        String assemblyBases1 = refBases1 + extBases1;

        // chr1:100:1 - links to chr2:50:-1
        JunctionAssembly assembly1 = createAssembly(CHR_1, 100, FORWARD, assemblyBases1, refBases1.length() - 1);

        // chr2:50:-1 - links to chr1:100:1 and faces chr2:200:1
        String refBases2 = refGenome.getBaseString(CHR_2, 50, 150);
        String extBases2 = refGenome.getBaseString(CHR_1, 51, 100); // was 101
        String assemblyBases2 = extBases2 + refBases2;

        JunctionAssembly assembly2 = createAssembly(CHR_2, 50, REVERSE, assemblyBases2, extBases2.length());

        // chr2:200:1 - links to chr3:200:-1 inverted to faces chr2:50
        String refBases3 = refGenome.getBaseString(CHR_2, 101, 200);
        String extBases3 = Nucleotides.reverseComplementBases(refGenome.getBaseString(CHR_3, 151, 201));
        String assemblyBases3 = refBases3 + extBases3;

        JunctionAssembly assembly3 = createAssembly(CHR_2, 200, FORWARD, assemblyBases3, refBases3.length() - 1);

        // chr3:200:1 - links to chr2:200:1 inverted and faces chr3:50:-1
        String refBases4 = refGenome.getBaseString(CHR_3, 101, 200);
        String extBases4 = refGenome.getBaseString(CHR_2, 151, 200);
        String assemblyBases4 = refBases4 + Nucleotides.reverseComplementBases(extBases4);

        JunctionAssembly assembly4 = createAssembly(CHR_3, 200, FORWARD, assemblyBases4, refBases4.length() - 1);

        // chr3:50:-1 - links to chr1:250:-1 inverted and faces chr3:200:1
        String refBases5 = refGenome.getBaseString(CHR_3, 50, 150);;
        String extBases5 = Nucleotides.reverseComplementBases(refGenome.getBaseString(CHR_1, 250, 300));
        String assemblyBases5 = extBases5 + refBases5;

        JunctionAssembly assembly5 = createAssembly(CHR_3, 50, REVERSE, assemblyBases5, extBases5.length());

        // chr1:250:-1 - links to chr3:50:-1
        String refBases6 = refGenome.getBaseString(CHR_1, 250, 350);
        String extBases6 = refGenome.getBaseString(CHR_3, 50, 100);
        String assemblyBases6 = Nucleotides.reverseComplementBases(extBases6) + refBases6;

        JunctionAssembly assembly6 = createAssembly(CHR_1, 250, REVERSE, assemblyBases6, extBases6.length());

        // now add junction, mate and discordant reads
        // each junction read has a mate either on the ref side or across a facing link to another segment
        String posJuncReadCigar = "100M50S";
        String negJuncReadCigar = "50S100M";
        String mateReadCigar = "100M";

        List<Read> candidates1 = Lists.newArrayList();
        List<Read> candidates2 = Lists.newArrayList();
        List<Read> candidates3 = Lists.newArrayList();
        List<Read> candidates4 = Lists.newArrayList();
        List<Read> candidates5 = Lists.newArrayList();
        List<Read> candidates6 = Lists.newArrayList();


        String juncBases = refGenome.getBaseString(CHR_1, 51, 150);

        Read juncRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 51, juncBases, posJuncReadCigar, CHR_2, 101, true);

        assembly1.addJunctionRead(juncRead1);

        // a mate rate in the next segment
        String mateBases = refGenome.getBaseString(CHR_2, 101, 200);
        Read mateRead = createRead(
                juncRead1.id(), CHR_2, 101, mateBases, mateReadCigar, CHR_1, 1, false);
        mateRead.bamRecord().setReadNegativeStrandFlag(true);
        setSecondInPair(mateRead.bamRecord());
        candidates2.add(mateRead);
        mateRead.setMateRead(juncRead1);
        juncRead1.setMateRead(mateRead);

        // a junction read and mate both in the first segment
        Read juncRead1b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 51, juncBases, posJuncReadCigar, CHR_1, 1, false);
        juncRead1b.bamRecord().setReadNegativeStrandFlag(true);

        assembly1.addJunctionRead(juncRead1b);

        // a mate rate in the same segment
        mateBases = refGenome.getBaseString(CHR_1, 1, 100);
        mateRead = createRead(
                juncRead1b.id(), CHR_1, 1, mateBases, mateReadCigar, CHR_1, 51, true);
        setSecondInPair(mateRead.bamRecord());
        candidates1.add(mateRead);
        mateRead.setMateRead(juncRead1b);
        juncRead1b.setMateRead(mateRead);

        // a junction read with its mate soft-clipping into the next junction
        Read juncRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 50, juncBases, negJuncReadCigar, CHR_2, 101, true);

        assembly2.addJunctionRead(juncRead2);

        mateBases = refGenome.getBaseString(CHR_2, 101, 200) + extBases3.substring(0, 10);
        mateRead = createRead(
                juncRead2.id(), CHR_2, 101, mateBases, "100M10S", CHR_2, 50, false);
        setSecondInPair(mateRead.bamRecord());
        mateRead.bamRecord().setReadNegativeStrandFlag(true);
        candidates2.add(mateRead);
        assembly3.addJunctionRead(mateRead);
        mateRead.setMateRead(juncRead2);
        juncRead2.setMateRead(mateRead);

        // an assembly 2 junction read with its mate in the following segment
        Read juncRead2b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 50, juncBases, negJuncReadCigar, CHR_3, 100, false);

        assembly2.addJunctionRead(juncRead2b);

        mateBases = refGenome.getBaseString(CHR_3, 100, 199);
        mateRead = createRead(
                juncRead2.id(), CHR_3, 100, mateBases, "100M", CHR_2, 50, false);
        setSecondInPair(mateRead.bamRecord());
        mateRead.bamRecord().setReadNegativeStrandFlag(true);
        candidates4.add(mateRead);
        mateRead.setMateRead(juncRead2b);
        juncRead2b.setMateRead(mateRead);

        // an assembly 3 junction read with its mate soft-clipped in the same section
        Read juncRead3 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 101, juncBases, posJuncReadCigar, CHR_2, 50, false);
        juncRead3.bamRecord().setReadNegativeStrandFlag(true);

        assembly3.addJunctionRead(juncRead3);

        mateBases = extBases2.substring(extBases2.length() - 10) + refGenome.getBaseString(CHR_2, 50, 149);
        mateRead = createRead(
                juncRead3.id(), CHR_3, 50, mateBases, "10S100M", CHR_2, 101, true);
        setSecondInPair(mateRead.bamRecord());
        assembly2.addJunctionRead(mateRead);
        candidates3.add(mateRead);
        candidates2.add(juncRead3);
        mateRead.setMateRead(juncRead3);
        juncRead3.setMateRead(mateRead);

        // an assembly 3 junction read with its mate in the next section
        Read juncRead3b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 101, juncBases, posJuncReadCigar, CHR_2, 101, false);

        assembly3.addJunctionRead(juncRead3b);

        mateBases = refGenome.getBaseString(CHR_3, 100, 199);
        mateRead = createRead(
                juncRead3b.id(), CHR_3, 100, mateBases, "100M", CHR_2, 101, false);
        setSecondInPair(mateRead.bamRecord());
        candidates4.add(mateRead);
        mateRead.setMateRead(juncRead3b);
        juncRead3b.setMateRead(mateRead);


        // an assembly 4 junction read with its mate soft-clipped in the same segment
        Read juncRead4 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_3, 101, juncBases, posJuncReadCigar, CHR_3, 50, true);

        assembly4.addJunctionRead(juncRead4);

        mateBases = extBases5.substring(extBases5.length() -10) + refGenome.getBaseString(CHR_3, 50, 149);
        mateRead = createRead(
                juncRead4.id(), CHR_3, 50, mateBases, "10S100M", CHR_3, 101, false);
        mateRead.bamRecord().setReadNegativeStrandFlag(true);
        setSecondInPair(mateRead.bamRecord());
        assembly5.addJunctionRead(mateRead);
        candidates4.add(mateRead);
        candidates5.add(juncRead4);
        mateRead.setMateRead(juncRead4);
        juncRead4.setMateRead(mateRead);

        // an assembly 4 junction read with its mate in the final segment
        Read juncRead4b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_3, 101, juncBases, posJuncReadCigar, CHR_1, 300, true);

        assembly4.addJunctionRead(juncRead4b);

        mateBases = refGenome.getBaseString(CHR_1, 300, 350);
        mateRead = createRead(
                juncRead4b.id(), CHR_1, 300, mateBases, "101M", CHR_3, 101, false);
        mateRead.bamRecord().setReadNegativeStrandFlag(true);
        setSecondInPair(mateRead.bamRecord());
        candidates6.add(mateRead);
        mateRead.setMateRead(juncRead4b);
        juncRead4b.setMateRead(mateRead);


        // an assembly 5 junction read with its mate in the final segment
        Read juncRead5 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_3, 50, juncBases, negJuncReadCigar, CHR_1, 300, true);

        assembly5.addJunctionRead(juncRead5);

        mateBases = refGenome.getBaseString(CHR_1, 300, 350);
        mateRead = createRead(
                juncRead5.id(), CHR_1, 300, mateBases, "101M", CHR_3, 101, false);
        setSecondInPair(mateRead.bamRecord());
        mateRead.bamRecord().setReadNegativeStrandFlag(true);
        candidates6.add(mateRead);
        mateRead.setMateRead(juncRead5);
        juncRead5.setMateRead(mateRead);

        // an assembly 6 with its mate in the previous segment
        Read juncRead6 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 250, juncBases, negJuncReadCigar, CHR_3, 100, false);

        assembly6.addJunctionRead(juncRead6);

        mateBases = refGenome.getBaseString(CHR_3, 100, 200);
        mateRead = createRead(
                juncRead6.id(), CHR_3, 100, mateBases, "101M", CHR_1, 250, false);
        setSecondInPair(mateRead.bamRecord());
        candidates5.add(mateRead);
        mateRead.setMateRead(juncRead6);
        juncRead6.setMateRead(mateRead);

        // add discordant reads as candidates to check they are assigned as support correctly
        String discCigar = "50M";
        Read discRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 20, REF_BASES_400.substring(20, 70), discCigar, CHR_3, 100, true);
        Read discRead2 = createRead(
                discRead1.id(), CHR_3, 100, REF_BASES_400.substring(100, 150), discCigar, CHR_1, 20, false);
        setSecondInPair(discRead2.bamRecord());

        assertTrue(isValidSupportCoordsVsJunction(discRead1, assembly1.junction().isForward(), assembly1.junction().Position));
        candidates1.add(discRead1);

        assertTrue(isValidSupportCoordsVsJunction(discRead2, assembly4.junction().isForward(), assembly4.junction().Position));
        candidates4.add(discRead2);

        assertFalse(isValidSupportCoordsVsJunction(discRead2, assembly5.junction().isForward(), assembly5.junction().Position));
        candidates5.add(discRead2);

        Read discRead3 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_2, 100, REF_BASES_400.substring(100, 150), discCigar, CHR_1, 300, true);
        Read discRead4 = createRead(
                discRead3.id(), CHR_1, 300, REF_BASES_400.substring(300, 350), discCigar, CHR_2, 100, false);
        discRead4.bamRecord().setReadNegativeStrandFlag(true);
        setSecondInPair(discRead4.bamRecord());

        assertFalse(isValidSupportCoordsVsJunction(discRead3, assembly2.junction().isForward(), assembly2.junction().Position));
        assertTrue(isValidSupportCoordsVsJunction(discRead3, assembly3.junction().isForward(), assembly3.junction().Position));
        candidates3.add(discRead3);

        assertFalse(isValidSupportCoordsVsJunction(discRead4, assembly1.junction().isForward(), assembly1.junction().Position));
        assertTrue(isValidSupportCoordsVsJunction(discRead4, assembly6.junction().isForward(), assembly6.junction().Position));
        candidates6.add(discRead4);

        RefBaseExtender refBaseExtender = new RefBaseExtender();
        refBaseExtender.findAssemblyCandidateExtensions(assembly1, candidates1);
        refBaseExtender.findAssemblyCandidateExtensions(assembly2, candidates2);
        refBaseExtender.findAssemblyCandidateExtensions(assembly3, candidates3);
        refBaseExtender.findAssemblyCandidateExtensions(assembly4, candidates4);
        refBaseExtender.findAssemblyCandidateExtensions(assembly5, candidates5);
        refBaseExtender.findAssemblyCandidateExtensions(assembly6, candidates6);

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
        assertEquals(5, assembly2.supportCount());
        assertEquals(6, assembly3.supportCount());
        assertEquals(5, assembly4.supportCount());
        assertEquals(3, assembly5.supportCount());
        assertEquals(3, assembly6.supportCount());
        assertEquals(1, getSupportTypeCount(assembly1, DISCORDANT));
        assertEquals(1, getSupportTypeCount(assembly2, DISCORDANT));
        assertEquals(1, getSupportTypeCount(assembly3, DISCORDANT));
        assertEquals(2, getSupportTypeCount(assembly4, DISCORDANT));
        assertEquals(0, getSupportTypeCount(assembly5, DISCORDANT));
        assertEquals(2, getSupportTypeCount(assembly6, DISCORDANT));

        // check final sequence and read full sequence indices
        String fullSequence = refGenome.getBaseString(CHR_1, 1, 100)
                + refGenome.getBaseString(CHR_2, 50, 200)
                + Nucleotides.reverseComplementBases(refGenome.getBaseString(CHR_3, 50, 200))
                + refGenome.getBaseString(CHR_1, 250, 350);

        AssemblyAlignment assemblyAlignment = new AssemblyAlignment(phaseSet);

        assertEquals(fullSequence, assemblyAlignment.fullSequence());

        // now check the full sequence read indices
        checkSupportReadSequenceIndices(assembly1, 0, 150, FORWARD);
        checkSupportReadSequenceIndices(assembly2, 50, 261, FORWARD);
        checkSupportReadSequenceIndices(assembly3, 50, 301, FORWARD);
        checkSupportReadSequenceIndices(assembly4, 151, 411, REVERSE);
        checkSupportReadSequenceIndices(assembly5, 251, 452, REVERSE);
        checkSupportReadSequenceIndices(assembly6, 352, 552, FORWARD);
    }

    protected static void checkSupportReadSequenceIndices(
            final JunctionAssembly assembly, int seqIndexStart, int seqIndexEnd, final Orientation orientation)
    {
        for(SupportRead read : assembly.support())
        {
            assertEquals(orientation, read.fullAssemblyOrientation());
            assertTrue(positionsWithin(read.fullAssemblyIndexStart(), read.fullAssemblyIndexEnd(), seqIndexStart, seqIndexEnd));
        }
    }

    @Test
    public void testChr363Chain()
    {
        // similar to chromosome 3-6-3 in COLO - links are:
        // chr1:200:1 -> chr1:300:1 - local INV
        // chr1:100:-1 -> chr1:200:1 - short facing TI
        // chr1:100:-1 -> chr1:300:1 - longer facing TI
        // chr1:100:-1 -> chr2:100:-1 - BND via a remote ref match

        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_400);
        refGenome.RefGenomeMap.put(CHR_2, REF_BASES_600);

        // chr1:100:-1 - has reads soft-clipping at 200 and others going out further and soft-clipping at 300
        String refBases1 = refGenome.getBaseString(CHR_1, 100, 150);
        String extBases1 = refGenome.getBaseString(CHR_2, 100, 200);
        String assemblyBases1 = Nucleotides.reverseComplementBases(extBases1) + refBases1;

        JunctionAssembly assembly1 = createAssembly(CHR_1, 100, REVERSE, assemblyBases1, extBases1.length());

        // chr1:200:1
        String refBases2 = refGenome.getBaseString(CHR_1, 100, 200);
        String extBases2 = refGenome.getBaseString(CHR_1, 251, 301);
        String assemblyBases2 = refBases2 + Nucleotides.reverseComplementBases(extBases2);

        JunctionAssembly assembly2 = createAssembly(CHR_1, 200, FORWARD, assemblyBases2, refBases2.length() - 1);

        // chr1:300:1
        String refBases3 = refGenome.getBaseString(CHR_1, 251, 300);
        String extBases3 = refGenome.getBaseString(CHR_1, 151, 201);
        String assemblyBases3 = refBases3 + Nucleotides.reverseComplementBases(extBases3);

        JunctionAssembly assembly3 = createAssembly(CHR_1, 300, FORWARD, assemblyBases3, refBases3.length() - 1);

        List<Read> candidates1 = Lists.newArrayList();
        List<Read> candidates2 = Lists.newArrayList();
        List<Read> candidates3 = Lists.newArrayList();

        // first junction reads for each assembly
        Read juncRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, assemblyBases1, "50S50M", CHR_2, 150, true);
        juncRead1.bamRecord().setReadNegativeStrandFlag(true);
        juncRead1.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, "100M");

        String juncRead1Bases = extBases1.substring(51, 101) + refBases1;
        Read juncRead1b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, juncRead1Bases, "50S50M", CHR_1, 251, true);

        Read juncRead1c = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, juncRead1Bases, "50S50M", CHR_1, 151, true);
        juncRead1c.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, "50M50S");

        Read juncRead1d = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, juncRead1Bases, "50S50M", CHR_1, 150, true);

        assembly1.addJunctionRead(juncRead1d);

        String mateBases1 = refGenome.getBaseString(CHR_1, 150, 250);
        Read mate1d = createRead(
                juncRead1d.id(), CHR_1, 150, mateBases1, "100M", CHR_1, 100, false);
        mate1d.bamRecord().setReadNegativeStrandFlag(true);

        juncRead1d.setMateRead(mate1d);
        mate1d.setMateRead(juncRead1d);
        setSecondInPair(mate1d.bamRecord());

        assembly1.addJunctionRead(juncRead1d);
        candidates1.add(mate1d);

        // build sufficient support to warrant a branch
        for(int i = 0; i < 4; ++i)
        {
            Read juncRead = cloneRead(juncRead1d, READ_ID_GENERATOR.nextId());
            Read mateRead = cloneRead(mate1d, READ_ID_GENERATOR.nextId());

            juncRead.setMateRead(mateRead);
            mateRead.setMateRead(juncRead);

            assembly1.addJunctionRead(juncRead);
            candidates1.add(mateRead);
        }

        assembly1.addJunctionRead(juncRead1);
        assembly1.addJunctionRead(juncRead1b);
        assembly1.addJunctionRead(juncRead1c);

        candidates2.add(juncRead1);

        candidates3.add(juncRead1);
        candidates3.add(juncRead1b);

        juncRead1 = cloneRead(juncRead1, READ_ID_GENERATOR.nextId());
        assembly1.addJunctionRead(juncRead1);
        candidates2.add(juncRead1);
        candidates3.add(juncRead1);

        Read juncRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 151, assemblyBases2, "50M50S", CHR_2, 150, true);

        assembly2.addJunctionRead(juncRead2);

        Read juncRead2b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 151, assemblyBases2, "50M50S", CHR_1, 220, false);

        assembly2.addJunctionRead(juncRead2b);

        Read juncRead2c = createRead(
                juncRead1c.id(), CHR_1, 151, assemblyBases2, "50M50S", CHR_1, 100, false);
        juncRead2c.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, "50S50M");
        juncRead2c.bamRecord().setReadNegativeStrandFlag(true);

        assembly2.addJunctionRead(juncRead2c);

        juncRead1c.setMateRead(juncRead2c);
        juncRead2c.setMateRead(juncRead1c);
        setSecondInPair(juncRead2c.bamRecord());
        candidates1.add(juncRead2c);
        candidates2.add(juncRead1c);

        for(int i = 0; i < 4; ++i)
        {
            Read juncRead = cloneRead(juncRead1c, READ_ID_GENERATOR.nextId());
            Read mateRead = cloneRead(juncRead2c, READ_ID_GENERATOR.nextId());

            assembly1.addJunctionRead(juncRead);
            assembly2.addJunctionRead(mateRead);

            juncRead.setMateRead(mateRead);
            mateRead.setMateRead(juncRead);
            candidates1.add(mateRead);
            candidates2.add(juncRead);
        }

        Read juncRead3 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 251, assemblyBases3, "50M50S", CHR_2, 150, true);

        assembly3.addJunctionRead(juncRead3);

        // the mate of assembly junc read 1b
        Read juncRead3b = createRead(
                juncRead1b.id(), CHR_1, 251, assemblyBases3, "50M50S", CHR_1, 100, true);

        juncRead3b.setMateRead(juncRead1b);
        juncRead1b.setMateRead(juncRead3b);
        setSecondInPair(juncRead3b.bamRecord());

        assembly3.addJunctionRead(juncRead3b);
        candidates1.add(juncRead3b);

        // a junction read spanning both ends of the short TI - counts towards both assembly 1 and 2
        String readBases = Nucleotides.reverseComplementBases(refGenome.getBaseString(CHR_2, 100, 110))
                + REF_BASES_400.substring(100, 201) + REF_BASES_400.substring(300, 310);
        Read sharedRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, readBases, "10S101M10S", CHR_2, 150, true);
        sharedRead1.bamRecord().setReadNegativeStrandFlag(true);

        assembly1.addJunctionRead(sharedRead1);
        assembly2.addJunctionRead(sharedRead1);

        // a junction mate covering the short TI junction
        readBases = REF_BASES_400.substring(151, 201) + REF_BASES_400.substring(251, 301);
        Read sharedRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 151, readBases, "50M50S", CHR_2, 150, true);

        candidates1.add(sharedRead2);
        assembly2.addJunctionRead(sharedRead2);

        // a junction mate extending past the short TI junction
        readBases = REF_BASES_400.substring(151, 251);
        Read mate1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 151, readBases, "100M", CHR_2, 150, true);
        candidates1.add(mate1);

        // discordant read hitting the longer TI junction
        readBases = REF_BASES_400.substring(251, 301) + Nucleotides.reverseComplementBases(REF_BASES_400.substring(191, 201));
        Read sharedRead3 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 251, readBases, "50M10S", CHR_2, 150, true);

        candidates1.add(sharedRead3);
        assembly3.addJunctionRead(sharedRead3);

        // junction mate linked to assembly 2
        readBases = REF_BASES_400.substring(220, 270);
        Read mate2 = createRead(
                juncRead2b.id(), CHR_1, 220, readBases, "50M", CHR_1, 151, false);

        mate2.setMateRead(juncRead2b);
        juncRead2b.setMateRead(mate2);
        setSecondInPair(mate2.bamRecord());

        candidates2.add(mate2);
        candidates3.add(mate2);

        RefBaseExtender refBaseExtender = new RefBaseExtender();
        refBaseExtender.findAssemblyCandidateExtensions(assembly1, candidates1);
        refBaseExtender.findAssemblyCandidateExtensions(assembly2, candidates2);
        refBaseExtender.findAssemblyCandidateExtensions(assembly3, candidates3);

        // set up the details for the remote ref assembly on chr 2
        RemoteRegionAssembler remoteRegionAssembler = new RemoteRegionAssembler(refGenome, null);

        Read remoteRead1 = createRead(
                juncRead1.id(), CHR_2, 150, refGenome.getBaseString(CHR_2, 150, 250), "100M",
                CHR_1, juncRead1.alignmentStart(), juncRead1.negativeStrand());

        Read remoteRead2 = createRead(sharedRead1.id(), CHR_2, 150, refGenome.getBaseString(CHR_2, 150, 250), "100M",
                CHR_1, sharedRead1.alignmentStart(), sharedRead1.negativeStrand());

        RemoteRegion remoteRegion = new RemoteRegion(
                new ChrBaseRegion(CHR_2, remoteRead1.alignmentStart(), remoteRead1.alignmentEnd()), remoteRead1.id(), RemoteReadType.DISCORDANT);

        remoteRegion.addReadDetails(remoteRead2.id(), remoteRead2.alignmentStart(), remoteRead2.alignmentEnd(), RemoteReadType.DISCORDANT);

        PhaseGroup phaseGroup = new PhaseGroup(assembly1, assembly2);
        phaseGroup.addAssembly(assembly3);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(refGenome, remoteRegionAssembler, phaseGroup);
        phaseSetBuilder.buildPhaseSets();

        assertEquals(1, phaseGroup.phaseSets().size());
        PhaseSet phaseSet = phaseGroup.phaseSets().get(0);

        assertEquals(3, phaseSet.assemblyLinks().size());
        assertEquals(4, phaseSet.assemblies().size());

        JunctionAssembly branchedAssembly = phaseSet.assemblies().stream()
                .filter(x -> x != assembly1 && x != assembly2 && x != assembly3)
                .findFirst().orElse(null);

        assertNotNull(branchedAssembly);

        assertEquals(201, assembly1.refBaseLength());
        assertEquals(201, assembly3.refBaseLength());
        assertEquals(101, assembly2.refBaseLength());
        assertEquals(101, branchedAssembly.refBaseLength());

        hasAssemblyLink(phaseSet.assemblyLinks(), assembly2, assembly3, LinkType.SPLIT);
        hasAssemblyLink(phaseSet.assemblyLinks(), assembly1, assembly3, LinkType.FACING);
        hasAssemblyLink(phaseSet.assemblyLinks(), branchedAssembly, assembly2, LinkType.FACING);

        assertEquals(20, assembly1.supportCount());
        assertEquals(14, assembly2.supportCount());
        assertEquals(5, assembly3.supportCount());
        assertEquals(0, getSupportTypeCount(assembly1, DISCORDANT));
        assertEquals(0, getSupportTypeCount(assembly2, DISCORDANT));

        String bases1 = Nucleotides.reverseComplementBases(refGenome.getBaseString(CHR_2, 100, 200));
        String bases2 = refGenome.getBaseString(CHR_1, 100, 200);
        String bases3 = Nucleotides.reverseComplementBases(refGenome.getBaseString(CHR_1, 100, 300));

        // the bases from the repeated link are now omitted
        String expectedFullSequence = bases1 + bases2 + "C" + bases3;

        AssemblyAlignment assemblyAlignment = new AssemblyAlignment(phaseSet);
        String actualFullSequence = assemblyAlignment.fullSequence();

        assertEquals(expectedFullSequence, actualFullSequence);

        /*
        checkSupportReadSequenceIndices(assembly1, 0, 150, FORWARD);
        checkSupportReadSequenceIndices(assembly2, 0, 150, FORWARD);
        checkSupportReadSequenceIndices(assembly3, 0, 150, FORWARD);
        */

    }
}
