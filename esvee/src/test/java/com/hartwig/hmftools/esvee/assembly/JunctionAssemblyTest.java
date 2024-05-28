package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID_2;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.recordSoftClipsAndCrossesJunction;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import org.junit.Test;

public class JunctionAssemblyTest
{
    @Test
    public void testMiscAssemblyFunctions()
    {
        assertEquals(0, mismatchesPerComparisonLength(1));
        assertEquals(1, mismatchesPerComparisonLength(10));
        assertEquals(2, mismatchesPerComparisonLength(100));
    }

    @Test
    public void testPosJunctionExtensionSequence()
    {
        String refBases = REF_BASES_200.substring(0, 20);
        String extBases = REF_BASES_200.substring(100, 140);

        Junction junction = new Junction(CHR_1, 29, FORWARD);

        // first a basic assembly with all reads agreeing
        String readBases = refBases + extBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        readBases = refBases + extBases.substring(0, 35);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        readBases = refBases + extBases.substring(0, 32);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        // a read without soft-clip but mismatches agreeing with the extension
        String alignedMatchingBases = readBases.substring(0, 22);
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), 10, alignedMatchingBases, makeCigarString(alignedMatchingBases, 0, 0));
        read4.bamRecord().setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        assertTrue(recordSoftClipsAndCrossesJunction(read4, junction));

        // similar but too long and matching the ref
        readBases = REF_BASES_200.substring(0, 30);
        Read read5 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, 0));
        assertFalse(recordSoftClipsAndCrossesJunction(read5, junction));

        List<Read> reads = List.of(read1, read2, read3, read4);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());

        assertEquals(36, extSeqBuilder.minSupportLength());

        String sequence = refBases.substring(19) + extBases;
        assertEquals(sequence, extSeqBuilder.junctionSequence());

        assertEquals(4, extSeqBuilder.formAssemblySupport().size());

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
        read4 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        reads = List.of(read1, read2, read3, read4);
        extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

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
        Junction junction = new Junction(CHR_1, juncPosition, REVERSE);

        // first a basic assembly with all reads agreeing

        String readBases = extBases + refBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), juncPosition, readBases, makeCigarString(readBases, extBases.length(), 0));

        readBases = extBases.substring(5) + refBases;
        int scLength = extBases.length() - 5;
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), juncPosition, readBases, makeCigarString(readBases, scLength, 0));

        readBases = extBases.substring(8) + refBases;
        scLength = extBases.length() - 8;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), juncPosition, readBases, makeCigarString(readBases, scLength, 0));

        String alignedMatchingBases = readBases.substring(readBases.length() - 22);
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), juncPosition - 2, alignedMatchingBases, makeCigarString(alignedMatchingBases, 0, 0));
        read4.bamRecord().setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        assertTrue(recordSoftClipsAndCrossesJunction(read4, junction));

        List<Read> reads = List.of(read1, read2, read3, read4);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());

        assertEquals(36, extSeqBuilder.minSupportLength());

        String sequence = extBases + refBases.substring(0, 1);
        assertEquals(sequence, extSeqBuilder.junctionSequence());

        assertEquals(4, extSeqBuilder.formAssemblySupport().size());

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
        read4 = createRead(READ_ID_GENERATOR.nextId(), juncPosition, readBases, makeCigarString(readBases, extBases.length(), 0));

        reads = List.of(read1, read2, read3, read4);
        extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());
        assertEquals(sequence, extSeqBuilder.junctionSequence());
        assertEquals(3, extSeqBuilder.formAssemblySupport().size());
    }

    @Test
    public void tesHighRepeatExtensionSequence()
    {
        String refBases = REF_BASES_200.substring(0, 20);

        Junction junction = new Junction(CHR_1, 29, FORWARD);

        String consensusExtBases = "AAAAAAAACCGTGTGTCCCAGTAGTAGTCCTTTT";
        String readBases = refBases + consensusExtBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, consensusExtBases.length()));

        // start with 2 reads agreeing to establish a base-line for repeats
        Read read1b = cloneRead(read1, READ_ID_GENERATOR.nextId());
        Read read1c = cloneRead(read1, READ_ID_GENERATOR.nextId());
        Read read1d = cloneRead(read1, READ_ID_GENERATOR.nextId());

        String extBases = "AAAAAAAAACCGTGTGTCCCAGTAGTAGTCCTTTT"; // extra A
        readBases = refBases + extBases;
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        extBases = "AAAAAACCGTGTGTCCCAGTAGTAGTCCTTTT"; // less As
        readBases = refBases + extBases;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        extBases = "AAAAAAAACCGTGTGTCCCAGTAGTCCTTTT"; // less AGT
        readBases = refBases + extBases;
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        extBases = "AAAAAAAACCGTGTGTGTGTCCCAGTAGTAGTCCTTTT"; // extra GTs
        readBases = refBases + extBases;
        Read read5 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        extBases = "AAAAAAAAACCGTGTGTGTCCCAGTAGTCCTTTT"; // 2+ mismatches
        readBases = refBases + extBases;
        Read read6 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        List<Read> reads = List.of(read1, read1b, read1c, read1d, read2, read3, read4, read5, read6);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());

        assertEquals(consensusExtBases.length() + 1, extSeqBuilder.minSupportLength());

        String consensusSequence = extSeqBuilder.junctionSequence();
        String sequence = refBases.substring(19) + consensusExtBases;
        assertEquals(sequence, consensusSequence);

        List<SupportRead> supportReads = extSeqBuilder.formAssemblySupport();
        assertEquals(4, supportReads.size());
        assertEquals(4, supportReads.stream().filter(x -> x.mismatchCount() == 0).count());
    }

    @Test
    public void testExactAssemblyDedup()
    {
        String assemblyBases = REF_BASES_200.substring(0, 50);

        Junction posJunction = new Junction(CHR_1, 60, FORWARD);

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

        assemblies.forEach(x -> x.addJunctionRead(sharedRead));
        assemblies.forEach(x -> x.addJunctionRead(sharedRead2));

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

    @Test
    public void testJunctionAssemblySplitting()
    {
        String refSequence = REF_BASES_400;

        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refSequence);

        Junction junction = new Junction(CHR_1, 100, FORWARD);

        String juncReadBases1 = refSequence.substring(51, 101) + refSequence.substring(200, 250);

        Read juncRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 51, juncReadBases1, "50M50S", CHR_1, 1, false);

        List<Read> junctionReads = Lists.newArrayList(juncRead1);

        for(int i = 0; i < 9; ++i)
        {
            junctionReads.add(cloneRead(juncRead1, READ_ID_GENERATOR.nextId()));
        }

        String juncReadBases2 = refSequence.substring(51, 101) + refSequence.substring(270, 320);

        Read juncRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 51, juncReadBases2, "50M50S", CHR_1, 1, false);

        junctionReads.add(juncRead2);

        for(int i = 0; i < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT - 1; ++i)
        {
            junctionReads.add(cloneRead(juncRead2, READ_ID_GENERATOR.nextId()));
        }

        JunctionAssembler junctionAssembler = new JunctionAssembler(junction);
        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(junctionReads);

        assertEquals(2, assemblies.size());
        JunctionAssembly assembly = assemblies.get(0);
        assertEquals(10, assembly.supportCount());

        assembly = assemblies.get(1);
        assertEquals(5, assembly.supportCount());
    }
}
