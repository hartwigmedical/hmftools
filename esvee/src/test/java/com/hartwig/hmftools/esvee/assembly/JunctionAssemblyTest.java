package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID_2;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.expandReferenceBases;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

public class JunctionAssemblyTest
{
    @Test
    public void testPosJunctionExtensionSequence()
    {
        String refBases = REF_BASES_200.substring(0, 20);
        String extBases = REF_BASES_200.substring(100, 140);

        Junction junction = new Junction(CHR_1, 29, POS_ORIENT);

        // first a basic assembly with all reads agreeing
        String readBases = refBases + extBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        readBases = refBases + extBases.substring(0, 35);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        readBases = refBases + extBases.substring(0, 32);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        List<Read> reads = List.of(read1, read2, read3);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads, 2);

        assertTrue(extSeqBuilder.isValid());

        assertEquals(36, extSeqBuilder.minSupportLength());

        String sequence = refBases.substring(19) + extBases;
        assertEquals(sequence, extSeqBuilder.junctionSequence());

        assertEquals(3, extSeqBuilder.formAssemblySupport().size());

        // now test with low qual mismatches in all 3 three reads but still overall agreement
        read2.getBases()[30] = MockRefGenome.getNextBase(read1.getBases()[30]);
        read2.getBases()[31] = MockRefGenome.getNextBase(read1.getBases()[31]);
        read2.getBaseQuality()[30] = 20;
        read2.getBaseQuality()[31] = 20;

        read3.getBases()[32] = MockRefGenome.getNextBase(read1.getBases()[32]);
        read3.getBases()[33] = MockRefGenome.getNextBase(read1.getBases()[33]);
        read3.getBaseQuality()[32] = 20;
        read3.getBaseQuality()[33] = 20;

        readBases = refBases + extBases.substring(0, 29) + "TTTT" + extBases.substring(33, 38);
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        reads = List.of(read1, read2, read3, read4);
        extSeqBuilder = new ExtensionSeqBuilder(junction, reads, 2);

        assertTrue(extSeqBuilder.isValid());
        assertEquals(sequence, extSeqBuilder.junctionSequence());
        assertEquals(3, extSeqBuilder.formAssemblySupport().size());
    }

    @Test
    public void testNegJunctionExtensionSequence()
    {
        String refBases = REF_BASES_200.substring(0, 20);
        String extBases = REF_BASES_200.substring(100, 140);

        int juncPosition = 100;
        Junction junction = new Junction(CHR_1, juncPosition, NEG_ORIENT);

        // first a basic assembly with all reads agreeing

        String readBases = extBases + refBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), juncPosition, readBases, makeCigarString(readBases, extBases.length(), 0));

        readBases = extBases.substring(5) + refBases;
        int scLength = extBases.length() - 5;
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), juncPosition, readBases, makeCigarString(readBases, scLength, 0));

        readBases = extBases.substring(8) + refBases;
        scLength = extBases.length() - 8;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), juncPosition, readBases, makeCigarString(readBases, scLength, 0));

        List<Read> reads = List.of(read1, read2, read3);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads, 2);

        assertTrue(extSeqBuilder.isValid());

        assertEquals(36, extSeqBuilder.minSupportLength());

        String sequence = extBases + refBases.substring(0, 1);
        assertEquals(sequence, extSeqBuilder.junctionSequence());

        assertEquals(3, extSeqBuilder.formAssemblySupport().size());

        // low qual mismatches in all 3 three reads but still overall agreement
        read2.getBases()[10] = MockRefGenome.getNextBase(read1.getBases()[30]);
        read2.getBases()[11] = MockRefGenome.getNextBase(read1.getBases()[31]);
        read2.getBaseQuality()[10] = 20;
        read2.getBaseQuality()[11] = 20;

        read3.getBases()[12] = MockRefGenome.getNextBase(read1.getBases()[32]);
        read3.getBases()[13] = MockRefGenome.getNextBase(read1.getBases()[33]);
        read3.getBaseQuality()[12] = 20;
        read3.getBaseQuality()[13] = 20;

        readBases = extBases.substring(0, 10) + "TTTT" + extBases.substring(14) + refBases;
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), juncPosition, readBases, makeCigarString(readBases, extBases.length(), 0));

        reads = List.of(read1, read2, read3, read4);
        extSeqBuilder = new ExtensionSeqBuilder(junction, reads, 2);

        assertTrue(extSeqBuilder.isValid());
        assertEquals(sequence, extSeqBuilder.junctionSequence());
        assertEquals(3, extSeqBuilder.formAssemblySupport().size());
    }

    @Test
    public void tesHighRepeatExtensionSequence()
    {
        String refBases = REF_BASES_200.substring(0, 20);

        Junction junction = new Junction(CHR_1, 29, POS_ORIENT);

        String consensusExtBases = "AAAAAAAACCGTGTGTCCAGTAGTAGTCCTTTT";
        String readBases = refBases + consensusExtBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, consensusExtBases.length()));

        // start with 2 reads agreeing to establish a base-line for repeats
        Read read1b = cloneRead(read1, READ_ID_GENERATOR.nextId());

        String extBases = "AAAAAAAAACCGTGTGTCCAGTAGTCCTTTT"; // extra A, less AGT
        readBases = refBases + extBases;
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        extBases = "AAAAAACCGTGTGTGTGTCCAGTAGTAGTCCTTTT"; // less As, extra GTs
        readBases = refBases + extBases;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        List<Read> reads = List.of(read1, read1b, read2, read3);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads, 2);

        assertTrue(extSeqBuilder.isValid());

        assertEquals(consensusExtBases.length() + 1, extSeqBuilder.minSupportLength());

        String sequence = refBases.substring(19) + consensusExtBases + "TT";
        assertEquals(sequence, extSeqBuilder.junctionSequence());

        assertEquals(4, extSeqBuilder.formAssemblySupport().size());
        assertEquals(2, extSeqBuilder.formAssemblySupport().stream().filter(x -> x.mismatchCount() == 2).count());
    }

    @Test
    public void testBuildJunctionSequencePosOrientation()
    {
        String refBases = REF_BASES_RANDOM_100.substring(0, 20) + "TTGGCCAATT";

        //   pos           10        20         30
        //   pos 012345678901234567890123456789 0123456789
        //       GATCGATCGATCGATCGATCTTGGCCAATT AACCGGGG (first sequence)
        // 1st-asm index  012345678901234567890 123456

        //       GATCGATCGATCGATCGATCTTGGCCAATT AATCGGTT (2nd sequence)
        // 2nd-asm index  012345678901234567890 12345678

        Junction junction = new Junction(CHR_1, 29, POS_ORIENT);

        // read 1 defines the first sequence
        Read read1 = createRead("READ_01", 10, refBases.substring(10, 30) + "AACCGGGG", "20M8S");
        Read read1b = cloneRead(read1, read1.getName() + "b");

        // read 2 has one mismatch on the first sequence and has a ref base mismatch
        String readRefBases = REF_BASES_RANDOM_100.substring(0, 20) + "TGGGCCAATT";
        Read read2 = createRead("READ_02", 9, readRefBases.substring(9, 30) + "ACCCGG", "21M6S");

        // read 3 defines the second sequence
        Read read3 = createRead("READ_03", 15, refBases.substring(15, 30) + "AATCGGTT", "15M8S");

        // read 4 matches has mismatches against both sequences, but only one with the second
        Read read4 = createRead("READ_04", 10, refBases.substring(10, 30) + "AATCGGA", "20M7S");

        // read 5 matches the second sequence and also has a ref base mismatch
        readRefBases = REF_BASES_RANDOM_100.substring(0, 20) + "TTGCCCAATT";
        Read read5 = createRead("READ_05", 10, readRefBases.substring(10, 30) + "AATCGGT", "20M7S");

        // has 1 mismatch against each so will support both
        Read read6 = createRead("READ_06", 10, refBases.substring(10, 30) + "AATCGGG", "20M7S");

        JunctionAssembly junctionSequence = AssemblyUtils.buildFromJunctionReads(
                junction, List.of(read1, read1b, read2, read3, read4, read5, read6), true);

        assertEquals(0, junctionSequence.junctionIndex());
        assertEquals(29, junctionSequence.minAlignedPosition());
        assertEquals(37, junctionSequence.maxAlignedPosition());

        assertEquals(5, junctionSequence.mismatches().allBaseMismatches().size());
        assertEquals(7, junctionSequence.supportCount());
        assertEquals(4, junctionSequence.mismatches().positionCount());

        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(junctionSequence);
        List<JunctionAssembly> allSequences = splitter.splitOnMismatches(5, 2);
        assertEquals(2, allSequences.size());

        assertEquals(4, allSequences.get(0).supportCount());
        assertEquals(4, allSequences.get(1).supportCount());

        allSequences.forEach(x -> expandReferenceBases(x));

        // for now support count is duplicated across reads for ref bases and post-junction bases
        assertEquals(4, allSequences.get(0).supportCount());
        assertEquals(1, allSequences.get(0).mismatches().positionCount());
        assertEquals(4, allSequences.get(1).supportCount());
        assertEquals(1, allSequences.get(1).mismatches().positionCount());
    }

    @Test
    public void testBuildJunctionSequenceNegOrientation()
    {
        String refBases = REF_BASES_RANDOM_100.substring(0, 20) + "TTGGCCAATT" + REF_BASES_RANDOM_100.substring(30, 50);

        //   pos           10        20         30
        //   pos 012345678901234567890123456789 0123456789
        //                  AACCGGGG TTGGCCAATT GATCGATCGATCGATCGATCTTGGCCAATT (2nd sequence)
        // 1st-asm index   012345678 9012345678 90123456

        //               TTTTCCTTGG TTGGCCAATT GATCGATCGATCGATCGATCTTGGCCAATT (1st sequence, designated since longest / high-qual SC)
        // 2nd-asm index   01234567 8901234567890 12345678

        Junction junction = new Junction(CHR_1, 20, NEG_ORIENT);

        // read 1 defines the first sequence
        Read read1 = createRead("READ_01", 20, "AACCGGGG" + refBases.substring(20, 32), "8S12M");

        // read 2 supports the first sequence
        Read read2 = createRead("READ_02", 20, "ACCGGGGT" + refBases.substring(21, 40), "7S20M");

        // read 3 defines the second sequence
        Read read3 = createRead("READ_03", 20, "TTTTCCTTGG" + refBases.substring(20, 35), "10S15M");

        // read 4 matches has 1 mismatch against both sequences
        Read read4 = createRead("READ_04", 20, "CCTGGG" + refBases.substring(20, 40), "6S20M");

        // read 5 matches the second sequence
        Read read5 = createRead("READ_05", 20, "TTCCTTGG" + refBases.substring(20, 37), "8S17M");

        JunctionAssembly junctionSequence = AssemblyUtils.buildFromJunctionReads(
                junction, List.of(read1, read2, read3, read4, read5), true);

        assertEquals(10, junctionSequence.junctionIndex());
        assertEquals(10, junctionSequence.minAlignedPosition());
        assertEquals(20, junctionSequence.maxAlignedPosition());

        assertEquals(4, junctionSequence.mismatches().allBaseMismatches().size());
        assertEquals(5, junctionSequence.supportCount());
        assertEquals(4, junctionSequence.mismatches().positionCount());

        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(junctionSequence);
        List<JunctionAssembly> allSequences = splitter.splitOnMismatches(5, 2);
        assertEquals(2, allSequences.size());

        assertEquals(3, allSequences.get(0).supportCount());
        assertEquals(3, allSequences.get(1).supportCount());

        allSequences.forEach(x -> expandReferenceBases(x));

        // for now support count is duplicated across reads for ref bases and post-junction bases
        JunctionAssembly firstSequence = allSequences.stream().filter(x -> x.initialRead() == read3).findFirst().orElse(null);

        assertEquals(3, firstSequence.supportCount());
        assertEquals(10, firstSequence.minAlignedPosition());
        assertEquals(39, firstSequence.maxAlignedPosition());
        assertEquals(0, firstSequence.mismatches().positionCount());

        String firstBases = firstSequence.formJunctionSequence(5);
        assertEquals("TTTTCCTTGG" + refBases.substring(20, 25), firstBases);

        JunctionAssembly secondSequence = allSequences.stream().filter(x -> x.initialRead() == read1).findFirst().orElse(null);

        assertEquals(3, secondSequence.supportCount());
        assertEquals(12, secondSequence.minAlignedPosition());
        assertEquals(39, secondSequence.maxAlignedPosition());
        assertEquals(0, secondSequence.mismatches().positionCount());

        String secondBases = secondSequence.formJunctionSequence(5);
        assertEquals("AACCGGGG" + refBases.substring(20, 25), secondBases);
    }

    @Test
    public void testExactAssemblyDedup()
    {
        String assemblyBases = REF_BASES_200.substring(0, 50);

        Junction posJunction = new Junction(CHR_1, 60, POS_ORIENT);

        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());

        JunctionAssembly assembly1 = new JunctionAssembly(
                posJunction, assemblyBases.getBytes(), baseQuals, 21); // was 40 - 89

        JunctionAssembly assembly2 = new JunctionAssembly(
                posJunction, assemblyBases.getBytes(), baseQuals, 21);

        JunctionAssembly assembly3 = new JunctionAssembly(
                posJunction, assemblyBases.getBytes(), baseQuals, 21);

        JunctionAssembly assembly4 = new JunctionAssembly(
                posJunction, assemblyBases.getBytes(), baseQuals, 21);

        JunctionAssembly assembly5 = new JunctionAssembly(
                posJunction, assemblyBases.getBytes(), baseQuals, 21);

        List<JunctionAssembly> assemblies = Lists.newArrayList(assembly1, assembly2, assembly3, assembly4, assembly5);

        Read sharedRead = createRead(TEST_READ_ID, 40, assemblyBases.substring(0, 30), "30M");
        Read sharedRead2 = createRead(TEST_READ_ID_2, 40, assemblyBases.substring(0, 30), "30M");

        assemblies.forEach(x -> x.addJunctionRead(sharedRead, false));
        assemblies.forEach(x -> x.addJunctionRead(sharedRead2, false));

        // add 2 mismatches to each of the other assemblies
        assembly2.bases()[2] = getNextBase(assembly2.bases()[2]);
        assembly2.bases()[30] = getNextBase(assembly2.bases()[30]);

        assembly3.bases()[10] = getNextBase(assembly2.bases()[10]);
        assembly3.bases()[20] = getNextBase(assembly2.bases()[20]);

        // add 5+ to the next one
        int index = 20;
        assembly4.bases()[index] = getNextBase(assembly4.bases()[index]);
        assembly4.bases()[++index] = getNextBase(assembly4.bases()[index]);
        assembly4.bases()[++index] = getNextBase(assembly4.bases()[index]);
        assembly4.bases()[++index] = getNextBase(assembly4.bases()[index]);
        assembly4.bases()[++index] = getNextBase(assembly4.bases()[index]);

        AssemblyDeduper.dedupJunctionAssemblies(assemblies);

        assertEquals(2, assemblies.size());
        assertEquals(3, assemblies.get(0).mergedAssemblyCount());
    }
}
