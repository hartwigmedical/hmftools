package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_600;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.getSupportTypeCount;
import static com.hartwig.hmftools.esvee.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.esvee.assembly.AssemblyLinksTest.checkSupportReadSequenceIndices;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.hasAssemblyLink;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.EXTENSION;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION_MATE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.utils.Strings;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteRegionAssembler;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
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
    public void testUnamppedSglPair()
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

        PhaseGroup phaseGroup = new PhaseGroup(assembly1, assembly2);

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

        AssemblyAlignment assemblyAlignment = new AssemblyAlignment(0, phaseSet);

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

        Read juncRead1b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, assemblyBases1, "50S50M", CHR_1, 251, true);

        Read juncRead1c = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, assemblyBases1, "50S50M", CHR_1, 151, true);
        juncRead1c.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, "50M50S");

        Read juncRead1d = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, assemblyBases1, "50S50M", CHR_1, 150, true);

        assembly1.addJunctionRead(juncRead1d);

        Read mate1d = createRead(
                juncRead1d.id(), CHR_1, 150, assemblyBases1, "100M", CHR_1, 100, false);
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

        remoteRegionAssembler.addMatchedReads(List.of(remoteRead1, remoteRead2), remoteRegion);

        PhaseGroup phaseGroup = new PhaseGroup(assembly1, assembly2);
        phaseGroup.addAssembly(assembly3);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(refGenome, remoteRegionAssembler, phaseGroup);
        phaseSetBuilder.buildPhaseSets();

        assertEquals(1, phaseGroup.phaseSets().size());
        PhaseSet phaseSet = phaseGroup.phaseSets().get(0);

        assertEquals(4, phaseSet.assemblyLinks().size());
        assertEquals(5, phaseSet.assemblies().size());

        JunctionAssembly remoteAssembly = phaseSet.assemblies().get(3);
        JunctionAssembly branchedAssembly = phaseSet.assemblies().get(4);

        hasAssemblyLink(phaseSet.assemblyLinks(), assembly2, assembly3, LinkType.SPLIT);
        hasAssemblyLink(phaseSet.assemblyLinks(), assembly1, assembly3, LinkType.FACING);
        hasAssemblyLink(phaseSet.assemblyLinks(), branchedAssembly, assembly2, LinkType.FACING);

        /*
        assertEquals(3, assembly1.supportCount());
        assertEquals(2, assembly2.supportCount());
        assertEquals(3, assembly3.supportCount());
        assertEquals(1, getSupportTypeCount(assembly1, DISCORDANT));
        assertEquals(0, getSupportTypeCount(assembly2, DISCORDANT));
        */

        String fullSequence = Nucleotides.reverseComplementBases(refGenome.getBaseString(CHR_2, 100, 249))
                + refGenome.getBaseString(CHR_1, 100, 300) + "A"
                + Nucleotides.reverseComplementBases(refGenome.getBaseString(CHR_1, 100, 200));

        AssemblyAlignment assemblyAlignment = new AssemblyAlignment(0, phaseSet);

        assertEquals(fullSequence, assemblyAlignment.fullSequence());

        /*
        checkSupportReadSequenceIndices(assembly1, 0, 150, FORWARD);
        checkSupportReadSequenceIndices(assembly2, 0, 150, FORWARD);
        checkSupportReadSequenceIndices(assembly3, 0, 150, FORWARD);
        */

    }
}
