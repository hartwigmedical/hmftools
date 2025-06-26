package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID_2;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.AssemblyDeduper.dedupProximateAssemblies;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readSoftClipsAndCrossesJunction;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
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
        assertEquals(1, mismatchesPerComparisonLength(15));
        assertEquals(1, mismatchesPerComparisonLength(100));
        assertEquals(2, mismatchesPerComparisonLength(101));
        assertEquals(3, mismatchesPerComparisonLength(201));
        assertEquals(3, mismatchesPerComparisonLength(399));
        assertEquals(4, mismatchesPerComparisonLength(500));
    }

    @Test
    public void testBasicJunctionAssemblies()
    {
        String extBases = REF_BASES_400.substring(300, 350);

        Junction junction = new Junction(CHR_1, 100, FORWARD);

        // assembly has 2 high-qual junction reads
        String readBases = REF_BASES_400.substring(51, 101) + extBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 51, readBases, makeCigarString(readBases, 0, extBases.length()));

        String readBases2 = readBases.substring(1);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 52, readBases, makeCigarString(readBases2, 0, extBases.length()));

        // the 3rd read doesn't overlap the junction enough to count as a junction read
        readBases = REF_BASES_400.substring(21, 101) + extBases.substring(0, 2);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 21, readBases, makeCigarString(readBases, 0, extBases.length()));

        List<Read> reads = List.of(read1, read2, read3);

        JunctionAssembler junctionAssembler = new JunctionAssembler(junction);

        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());

        JunctionAssembly assembly = assemblies.get(0);

        assertEquals(50, assembly.refBaseLength());
        assertEquals(51, assembly.refBasePosition());
        assertEquals(2, assembly.supportCount()); // read counts as support, just doesn't extend the ref bases
    }

    @Test
    public void testPosJunctionExtensionSequence()
    {
        String refBases = REF_BASES_200.substring(10, 30);
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

        MockRefGenome refGenome = new MockRefGenome();
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_200);

        read4.bamRecord().setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        assertTrue(readSoftClipsAndCrossesJunction(read4, junction, refGenome));

        // similar but matches the ref
        readBases = REF_BASES_200.substring(10, 40);
        Read read5 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, 0));
        assertFalse(readSoftClipsAndCrossesJunction(read5, junction, refGenome));

        List<Read> reads = List.of(read1, read2, read3, read4);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());

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

        List<Read> reads = List.of(read1, read2, read3, read4);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());

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

        String bufferBases = "CAGTAGGAGCTTTA";
        String consensusRead = "AAAAAAAA";
        String consensusExtBases = bufferBases + consensusRead + bufferBases;
        String readBases = refBases + consensusExtBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, consensusExtBases.length()));

        // start with 2 reads agreeing to establish a base-line for repeats
        Read read1b = cloneRead(read1, READ_ID_GENERATOR.nextId());
        Read read1c = cloneRead(read1, READ_ID_GENERATOR.nextId());

        String extBases = bufferBases + consensusRead + "A" + bufferBases; // extra A
        readBases = refBases + extBases;
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));
        extBases = bufferBases + consensusRead.substring(3) + bufferBases; // less As
        readBases = refBases + extBases;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        extBases = bufferBases + consensusRead + "AAA" + bufferBases; // extra As
        readBases = refBases + extBases;
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), 10, readBases, makeCigarString(readBases, 0, extBases.length()));

        List<Read> reads = List.of(read1, read1b, read1c, read2, read3, read4);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());

        String consensusSequence = extSeqBuilder.junctionSequence();
        String sequence = refBases.substring(19) + consensusExtBases;
        assertEquals(sequence, consensusSequence);

        List<SupportRead> supportReads = extSeqBuilder.formAssemblySupport();
        assertEquals(4, supportReads.size());
        assertEquals(4, supportReads.stream().filter(x -> x.extensionBaseMismatches() == 0).count());
    }

    @Test
    public void tesRepeatMismatchesExtensionSequence()
    {
        String refBases = REF_BASES_200.substring(0, 20);

        int junctionPosition = 50;
        Junction junction = new Junction(CHR_1, 50, REVERSE);

        String buffer = "GATCGTAGGATC";
        String caRepeat = "CACACACA";

        String extBases1 = buffer + caRepeat + buffer;
        String readBases = extBases1 + refBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), junctionPosition, readBases, makeCigarString(readBases, extBases1.length(), 0));

        String extBases2 = buffer + caRepeat + "CA" + buffer; // 1 extra
        readBases = extBases2 + refBases;
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), junctionPosition, readBases, makeCigarString(readBases, extBases2.length(), 0));

        String extBases3 = buffer + caRepeat.substring(2) + buffer; // one less CA
        readBases = extBases3 + refBases;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), junctionPosition, readBases, makeCigarString(readBases, extBases3.length(), 0));

        String extBases4 = buffer + caRepeat + "CACA" + buffer; // 2 extra CAs
        readBases = extBases4 + refBases;
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), junctionPosition, readBases, makeCigarString(readBases, extBases4.length(), 0));

        Read read1b = cloneRead(read1, READ_ID_GENERATOR.nextId());

        List<Read> reads = List.of(read1, read2, read3, read4, read1b);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());

        String consensusSequence = extSeqBuilder.junctionSequence();

        String consensusExtBases = buffer + caRepeat + buffer; // matches read 1

        String sequence = consensusExtBases + refBases.substring(0, 1);
        assertEquals(sequence, consensusSequence);

        List<SupportRead> supportReads = extSeqBuilder.formAssemblySupport();
        assertEquals(4, supportReads.size());
        assertEquals(4, supportReads.stream().filter(x -> x.extensionBaseMismatches() == 0).count());
    }

    @Test
    public void tesRefBaseRepeatExpansionMismatchesExtensionSequence()
    {
        int junctionPosition = 50;
        Junction junction = new Junction(CHR_1, junctionPosition, REVERSE);

        String refBases = "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG" + REF_BASES_200.substring(0, 20); // has 10 repeats
        String extraExtBases = "GATCGTAGGATCGTAGGATCGTAGGATC"; // other non-repeated bases

        String extBases1 = extraExtBases + "CAGCAGCAGCAG"; // has 4 repeats
        String readBases = extBases1 + refBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), junctionPosition, readBases, makeCigarString(readBases, extBases1.length(), 0));

        Read read1b = cloneRead(read1, READ_ID_GENERATOR.nextId());

        String extBases3 = extBases1 + "CAGCAG"; // 2 extra CAGs
        readBases = extBases3 + refBases;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), junctionPosition, readBases, makeCigarString(readBases, extBases3.length(), 0));

        String extBases4 = extraExtBases + "CAGCAG"; // 2 less CAGs
        readBases = extBases4 + refBases;
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), junctionPosition, readBases, makeCigarString(readBases, extBases4.length(), 0));

        String extBases5 = "CAG"; // insufficient overlap bases
        readBases = extBases5 + refBases;
        Read read5 = createRead(READ_ID_GENERATOR.nextId(), junctionPosition, readBases, makeCigarString(readBases, extBases5.length(), 0));

        List<Read> reads = List.of(read1, read1b, read3, read4, read5);

        ExtensionSeqBuilder extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());
        assertNotNull(extSeqBuilder.maxRepeat());
        assertEquals(4, extSeqBuilder.maxRepeat().Count);
        assertEquals(10, extSeqBuilder.refBaseRepeatCount());

        String consensusSequence = extSeqBuilder.junctionSequence();

        String consensusExtBases = extBases1; // matches read 1

        String sequence = consensusExtBases + refBases.substring(0, 1);
        assertEquals(sequence, consensusSequence);

        List<SupportRead> supportReads = extSeqBuilder.formAssemblySupport();
        assertEquals(4, supportReads.size());
        assertEquals(4, supportReads.stream().filter(x -> x.extensionBaseMismatches() == 0).count());

        String buildInfo = extSeqBuilder.buildInformation();
        assertEquals("5;5;0;true;28:CAG:4:10;2:2:1:0", buildInfo);

        // test a forward orientation junction
        junction = new Junction(CHR_1, junctionPosition, FORWARD);

        refBases = REF_BASES_200.substring(0, 20) + "CACACACACACACACACACA"; // has 10 repeats

        extBases1 = "CACACACA" + extraExtBases; // has 4 repeats
        readBases = refBases + extBases1;
        int readStartPos = junctionPosition - refBases.length() + 1;
        read1 = createRead(READ_ID_GENERATOR.nextId(), readStartPos, readBases, makeCigarString(readBases, 0, extBases1.length()));

        read1b = cloneRead(read1, READ_ID_GENERATOR.nextId());

        extBases3 = "CACA" + extBases1; // 2 extra repeats
        readBases = refBases + extBases3;
        read3 = createRead(READ_ID_GENERATOR.nextId(), readStartPos, readBases, makeCigarString(readBases, 0, extBases3.length()));

        extBases4 = "CACA" + extraExtBases; // 2 less CAs
        readBases = refBases + extBases4;
        read4 = createRead(READ_ID_GENERATOR.nextId(), readStartPos, readBases, makeCigarString(readBases, 0, extBases4.length()));

        extBases5 = "CA"; // insufficient overlap bases
        readBases = refBases + extBases5;
        read5 = createRead(READ_ID_GENERATOR.nextId(), readStartPos, readBases, makeCigarString(readBases, 0, extBases5.length()));

        reads = List.of(read1, read1b, read3, read4, read5);

        extSeqBuilder = new ExtensionSeqBuilder(junction, reads);

        assertTrue(extSeqBuilder.isValid());
        assertNotNull(extSeqBuilder.maxRepeat());
        assertEquals(4, extSeqBuilder.maxRepeat().Count);
        assertEquals(10, extSeqBuilder.refBaseRepeatCount());

        consensusSequence = extSeqBuilder.junctionSequence();

        consensusExtBases = extBases1; // matches read 1

        sequence = refBases.substring(refBases.length() - 1) + consensusExtBases;
        assertEquals(sequence, consensusSequence);

        supportReads = extSeqBuilder.formAssemblySupport();
        assertEquals(4, supportReads.size());
        assertEquals(4, supportReads.stream().filter(x -> x.extensionBaseMismatches() == 0).count());
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

        assertTrue(SequenceCompare.matchedAssemblySequences(assembly1, assembly2));
        assertTrue(SequenceCompare.matchedAssemblySequences(assembly1, assembly3));
        assertFalse(SequenceCompare.matchedAssemblySequences(assembly1, assembly4));
        assertTrue(SequenceCompare.matchedAssemblySequences(assembly1, assembly5));
    }

    @Test
    public void testAssemblyNoDedupWithSharedReads()
    {
        Junction indelJunctionNeg = new Junction(CHR_1, 130, REVERSE, false, true, false);
        Junction indelJunctionPos = new Junction(CHR_1, 131, FORWARD, false, true, false);

        String leftExtBases = REF_BASES_400.substring(200, 240);
        String rightExtBases = REF_BASES_400.substring(260, 300);
        String leftRefBases = REF_BASES_400.substring(100, 131);
        String rightRefBases = REF_BASES_400.substring(131, 171);
        String insertBases =rightRefBases; // a duplication

        String indelReadBases = leftExtBases + leftRefBases + insertBases + rightRefBases + rightExtBases;

        // both pairs of junctions share a read which supports them both, so are not deduped
        Read indelRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100,
                indelReadBases,
                "40S31M40I40M40S", CHR_1, 300, true);

        JunctionAssembly indelAssemblyPos = new JunctionAssembly(
                indelJunctionPos, indelRead.getBases(), indelRead.getBaseQuality(),
                leftExtBases.length() + leftRefBases.length());

        indelAssemblyPos.addJunctionRead(indelRead);
        indelAssemblyPos.setIndelCoords(indelRead.indelCoords());

        JunctionAssembly indelAssemblyNeg = new JunctionAssembly(
                indelJunctionNeg, indelRead.getBases(), indelRead.getBaseQuality(),
                leftExtBases.length() + leftRefBases.length() + 1);

        indelAssemblyNeg.addJunctionRead(indelRead);
        indelAssemblyNeg.setIndelCoords(indelRead.indelCoords());

        Junction splitJunctionPos = new Junction(CHR_1, 170, FORWARD);
        Junction splitJunctionNeg = new Junction(CHR_1, 100, REVERSE);

        JunctionAssembly splitAssemblyNeg = new JunctionAssembly(
                splitJunctionNeg, indelRead.getBases(), indelRead.getBaseQuality(), leftExtBases.length());
        splitAssemblyNeg.addJunctionRead(indelRead);

        JunctionAssembly splitAssemblyPos = new JunctionAssembly(
                splitJunctionPos, indelRead.getBases(), indelRead.getBaseQuality(),
                indelRead.basesLength() - rightExtBases.length() - 1);
        splitAssemblyPos.addJunctionRead(indelRead);

        List<JunctionAssembly> existingAssemblies = Lists.newArrayList(splitAssemblyPos, splitAssemblyNeg);
        List<JunctionAssembly> newAssemblies = Lists.newArrayList(indelAssemblyPos, indelAssemblyNeg);
        List<JunctionAssembly> dedupIndels = Lists.newArrayList();

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(2, existingAssemblies.size());
        assertEquals(2, newAssemblies.size());

        // again without a shared read and they are deduped
        Read splitRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100,
                indelReadBases, "40S31M30M70S", CHR_1, 300, true);

        splitAssemblyNeg = new JunctionAssembly(
                splitJunctionNeg, splitRead.getBases(), splitRead.getBaseQuality(), leftExtBases.length());
        splitAssemblyNeg.addJunctionRead(splitRead);

        splitAssemblyPos = new JunctionAssembly(
                splitJunctionPos, splitRead.getBases(), splitRead.getBaseQuality(),
                splitRead.basesLength() - rightExtBases.length() - 1);
        splitAssemblyPos.addJunctionRead(splitRead);

        existingAssemblies = Lists.newArrayList(splitAssemblyPos, splitAssemblyNeg);
        newAssemblies = Lists.newArrayList(indelAssemblyPos, indelAssemblyNeg);

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(2, newAssemblies.size());
        assertTrue(existingAssemblies.isEmpty()); // the indel assemblies are kept
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

        Read juncRead1b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 52, juncReadBases1.substring(1), "49M50S", CHR_1,
                1, false);

        junctionReads.add(juncRead1b);

        String juncReadBases2 = refSequence.substring(51, 101) + refSequence.substring(270, 320);

        Read juncRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 51, juncReadBases2, "50M50S", CHR_1, 1, false);

        junctionReads.add(juncRead2);

        for(int i = 0; i < ASSEMBLY_SPLIT_MIN_READ_SUPPORT - 1; ++i)
        {
            junctionReads.add(cloneRead(juncRead2, READ_ID_GENERATOR.nextId()));
        }

        Read juncRead2b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 52, juncReadBases2.substring(1), "49M50S", CHR_1,
                1, false);

        junctionReads.add(juncRead2b);

        JunctionAssembler junctionAssembler = new JunctionAssembler(junction);
        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(junctionReads);

        assertEquals(2, assemblies.size());
        JunctionAssembly assembly = assemblies.get(0);
        assertEquals(11, assembly.supportCount());

        assembly = assemblies.get(1);
        assertEquals(6, assembly.supportCount());
    }

    @Test
    public void testJunctionAssemblyRefSideSoftClipReadAssignment()
    {
        // reads which cross the dominant ref-side soft clip are removed from junction support
        Junction junction = new Junction(CHR_1, 200, FORWARD);

        String extBasesShared = "AAGGCCTTAA";
        String extBases1 = extBasesShared + REF_BASES_400.substring(300, 340);
        String extBases2 = extBasesShared + REF_BASES_400.substring(350, 390);

        // reads which support the first assembly and which soft-clip to the left
        String otherSoftClipBases = REF_BASES_400.substring(0, 50);
        String juncReadBases1 = otherSoftClipBases + REF_BASES_400.substring(151, 201) + extBases1;

        Read juncRead1 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 151, juncReadBases1, "50S50M50S", CHR_2, 100, false);

        Read juncRead1b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 152, juncReadBases1.substring(1), "50S49M50S", CHR_2,
                100, false);

        List<Read> junctionReads = Lists.newArrayList(juncRead1, juncRead1b);

        for(int i = 0; i < 5; ++i)
        {
            Read juncRead1c = cloneRead(juncRead1, READ_ID_GENERATOR.nextId());
            junctionReads.add(juncRead1c);
        }

        String juncReadBases2 = REF_BASES_400.substring(101, 201) + extBases2;
        Read juncRead2 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 101, juncReadBases2, "100M50S", CHR_2,
                100, false);

        Read juncRead2b = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 102, juncReadBases2.substring(1), "99M50S", CHR_2,
                100, false);

        junctionReads.add(juncRead2);
        junctionReads.add(juncRead2b);

        for(int i = 0; i < 5; ++i)
        {
            Read juncRead2c = cloneRead(juncRead2, READ_ID_GENERATOR.nextId());
            junctionReads.add(juncRead2c);
        }

        // then some reads which could support either extension but which go past the soft-clip and will be removed from the first assembly
        String juncReadBases3 = REF_BASES_400.substring(100, 201) + extBasesShared;
        Read juncRead3 = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, juncReadBases3, "101M10S", CHR_2,
                100, false);

        junctionReads.add(juncRead3);

        JunctionAssembler junctionAssembler = new JunctionAssembler(junction);
        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(junctionReads);

        assertEquals(2, assemblies.size());

        JunctionAssembly assembly = assemblies.get(0);
        assertEquals(8, assembly.supportCount());
        assertEquals(100, assembly.refBasePosition());

        assembly = assemblies.get(1);
        assertEquals(7, assembly.supportCount());
        assertEquals(151, assembly.refBasePosition());
    }
}
