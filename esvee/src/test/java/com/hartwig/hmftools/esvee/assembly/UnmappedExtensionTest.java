package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_600;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.setSecondInPair;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssembly;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.esvee.assembly.phase.PhaseSetBuilder;
import com.hartwig.hmftools.esvee.assembly.phase.RemoteReadExtractor;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;

import org.junit.Test;

public class UnmappedExtensionTest
{
    @Test
    public void testUnmappedReadExtension()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_400);

        String refBases = refGenome.getBaseString(CHR_1, 1, 101);
        String extBases = REF_BASES_400.substring(200, 250);
        String assemblyBases = refBases + extBases;

        JunctionAssembly assembly = createAssembly(CHR_1, 100, FORWARD, assemblyBases, refBases.length() - 1);

        // create a set of unmapped reads which will consider for extension
        String unmappedReadBases1 = REF_BASES_400.substring(200, 300);
        Read unmappedRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases1, NO_CIGAR, CHR_1, 1, true);
        unmappedRead1.bamRecord().setReadUnmappedFlag(true);

        // next read will need its bases reversed, and will only be added once the first one has extended the extension bases
        String unmappedReadBases2 = Nucleotides.reverseComplementBases(REF_BASES_400.substring(250, 350));
        Read unmappedRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases2, NO_CIGAR, CHR_1, 1, false);
        unmappedRead2.bamRecord().setReadUnmappedFlag(true);

        // next 3 reads require testing for consensus bases, and the first has too many mismatches to be addd
        String unmappedReadBases3 = REF_BASES_400.substring(300, 400);
        String unmappedReadBases3a = unmappedReadBases3.substring(0, 60) + "AAAAA" + unmappedReadBases3.substring(65);

        Read unmappedRead3a = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases3a, NO_CIGAR, CHR_1, 1, true);
        unmappedRead3a.bamRecord().setReadUnmappedFlag(true);

        // minor mismatches for the others
        String unmappedReadBases3b = unmappedReadBases3.substring(0, 70)
                + MockRefGenome.getNextBase(unmappedReadBases3.substring(70, 71)) + unmappedReadBases3.substring(71);
        Read unmappedRead3b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases3b, NO_CIGAR, CHR_1, 1, true);
        unmappedRead3b.bamRecord().setReadUnmappedFlag(true);

        String unmappedReadBases3c = Nucleotides.reverseComplementBases(unmappedReadBases3);
        Read unmappedRead3c = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases3c, NO_CIGAR, CHR_1, 1, false);
        unmappedRead3c.bamRecord().setReadUnmappedFlag(true);

        assembly.unmappedReads().add(unmappedRead1);
        assembly.unmappedReads().add(unmappedRead2);
        assembly.unmappedReads().add(unmappedRead3a);
        assembly.unmappedReads().add(unmappedRead3b);
        assembly.unmappedReads().add(unmappedRead3c);

        PhaseGroup phaseGroup = new PhaseGroup(assembly, null);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(refGenome, new RemoteReadExtractor(null), phaseGroup);
        phaseSetBuilder.buildPhaseSets();

        assertEquals(200, assembly.extensionLength());
        assertEquals(4, assembly.supportCount());
        assertFalse(assembly.hasReadSupport(unmappedRead3a));
    }

    @Test
    public void testUnmappedReadExtensionNegJunction()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_400);

        String refBases = refGenome.getBaseString(CHR_1, 300, 400);
        String extBases = REF_BASES_400.substring(150, 200);
        String assemblyBases = extBases + refBases;

        JunctionAssembly assembly = createAssembly(CHR_1, 300, REVERSE, assemblyBases, extBases.length());

        // 4 reads with varying levels of mismatches
        String unmappedReadBases1 = REF_BASES_400.substring(100, 200);
        Read unmappedRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases1, NO_CIGAR, CHR_1, 1, true);
        unmappedRead1.bamRecord().setReadUnmappedFlag(true);

        // next read will need its bases reversed, and will only be added once the first one has extended the extension bases
        String unmappedReadBases2 = Nucleotides.reverseComplementBases(REF_BASES_400.substring(100, 200));
        Read unmappedRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases2, NO_CIGAR, CHR_1, 1, false);
        unmappedRead2.bamRecord().setReadUnmappedFlag(true);

        String unmappedReadBases3 = REF_BASES_400.substring(100, 200);
        String unmappedReadBases3a = unmappedReadBases3.substring(0, 30) + "AAAAA" + unmappedReadBases3.substring(35);

        Read unmappedRead3a = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases3a, NO_CIGAR, CHR_1, 1, true);
        unmappedRead3a.bamRecord().setReadUnmappedFlag(true);

        // minor mismatches for the others
        String unmappedReadBases3b = unmappedReadBases3.substring(0, 10)
                + MockRefGenome.getNextBase(unmappedReadBases3.substring(10, 11)) + unmappedReadBases3.substring(11);
        Read unmappedRead3b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases3b, NO_CIGAR, CHR_1, 1, true);
        unmappedRead3b.bamRecord().setReadUnmappedFlag(true);

        String unmappedReadBases3c = unmappedReadBases3.substring(0, 40)
                + MockRefGenome.getNextBase(unmappedReadBases3.substring(40, 41)) + unmappedReadBases3.substring(41);
        Read unmappedRead3c = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, unmappedReadBases3c, NO_CIGAR, CHR_1, 1, true);
        unmappedRead3c.bamRecord().setReadUnmappedFlag(true);

        assembly.unmappedReads().add(unmappedRead1);
        assembly.unmappedReads().add(unmappedRead2);
        assembly.unmappedReads().add(unmappedRead3a);
        assembly.unmappedReads().add(unmappedRead3b);
        assembly.unmappedReads().add(unmappedRead3c);

        PhaseGroup phaseGroup = new PhaseGroup(assembly, null);

        PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(refGenome, new RemoteReadExtractor(null), phaseGroup);
        phaseSetBuilder.buildPhaseSets();

        assertEquals(100, assembly.extensionLength());
        assertEquals(4, assembly.supportCount());
        assertFalse(assembly.hasReadSupport(unmappedRead3a));
    }

}
