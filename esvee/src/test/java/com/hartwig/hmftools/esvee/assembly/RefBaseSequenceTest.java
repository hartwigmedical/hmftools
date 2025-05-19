package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.calcIndelInferredUnclippedPositions;
import static com.hartwig.hmftools.esvee.assembly.RefBaseSeqBuilder.readRefBaseLength;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

import org.junit.Test;

public class RefBaseSequenceTest
{
    @Test
    public void testForwardReadParseState()
    {
        String readBases = REF_BASES_200.substring(1, 41);
        Read read = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "10S30M");
        int readJunctionIndex = 10; // soft-clip length
        boolean isForwardJunction = false;

        RefReadParseState readState = new RefReadParseState(
                isForwardJunction, read, readJunctionIndex, readRefBaseLength(read, readJunctionIndex, isForwardJunction));

        assertEquals(11, readState.refPosition());
        assertEquals(10, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(11), (char)readState.currentBase());

        for(int i = 0; i < 29; ++i)
        {
            readState.moveNext();
        }

        assertEquals(40, readState.refPosition());
        assertEquals(39, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(40), (char)readState.currentBase());
        assertFalse(readState.exhausted());

        readState.moveNext();
        assertTrue(readState.exhausted());

        // test again with indels
        read = createRead(READ_ID_GENERATOR.nextId(), 1, REF_BASES_200.substring(1, 36), "5M5D5M5I10M10S");

        readState = new RefReadParseState(
                isForwardJunction, read, 0, readRefBaseLength(read, 0, isForwardJunction));

        assertEquals(1, readState.refPosition());
        assertEquals(0, readState.readIndex());
        assertEquals(M, readState.operator());

        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(6, readState.refPosition());
        assertEquals(4, readState.readIndex());
        assertEquals(D, readState.operator());

        // move through the delete
        readState.moveNext();

        assertEquals(7, readState.refPosition());
        assertEquals(4, readState.readIndex());
        assertEquals(D, readState.operator());

        readState.moveNext();
        readState.moveNext();
        readState.moveNext();

        assertEquals(10, readState.refPosition());
        assertEquals(4, readState.readIndex());
        assertEquals(D, readState.operator());

        readState.moveNext();

        assertEquals(11, readState.refPosition());
        assertEquals(5, readState.readIndex());
        assertEquals(M, readState.operator());

        // move onto insert
        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(15, readState.refPosition());
        assertEquals(10, readState.readIndex());
        assertEquals(I, readState.operator());

        readState.moveNext();
        assertEquals(15, readState.refPosition());
        assertEquals(11, readState.readIndex());
        assertEquals(I, readState.operator());

        // through rest of insert
        readState.moveNext();
        readState.moveNext();
        readState.moveNext();
        readState.moveNext();

        assertEquals(16, readState.refPosition());
        assertEquals(15, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(16), (char)readState.currentBase());

        // stop at final soft-clipping
        for(int i = 0; i < 10; ++i)
        {
            readState.moveNext();
        }

        assertTrue(readState.exhausted());
    }

    @Test
    public void testReverseReadState()
    {
        String readBases = REF_BASES_200.substring(1, 41);
        Read read = createRead(READ_ID_GENERATOR.nextId(), 1, readBases, "30M10S");
        int readJunctionIndex = 29;
        boolean isForwardJunction = true;

        RefReadParseState readState = new RefReadParseState(
                isForwardJunction, read, readJunctionIndex, readRefBaseLength(read, readJunctionIndex, isForwardJunction));

        assertEquals(30, readState.refPosition());
        assertEquals(29, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(30), (char)readState.currentBase());

        for(int i = 0; i < 29; ++i)
        {
            readState.moveNext();
        }

        assertEquals(1, readState.refPosition());
        assertEquals(0, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(1), (char)readState.currentBase());
        assertFalse(readState.exhausted());

        readState.moveNext();
        assertTrue(readState.exhausted());

        // test again with indels
        read = createRead(READ_ID_GENERATOR.nextId(), 11, REF_BASES_200.substring(1, 36), "10S10M5D5M5I5M");
        // bases 31-35, index 30-34 = 5M
        // bases 30, index 25-29 = 5I
        // bases 26-30, index 20-24 = 5M
        // bases 21-25,  = 5D
        // bases 11-20, index 10-19 = 10M
        // bases 1-10, index 0-9 = 10S

        readJunctionIndex = read.basesLength() - 1;
        readState = new RefReadParseState(
                isForwardJunction, read, readJunctionIndex, readRefBaseLength(read, readJunctionIndex, isForwardJunction));

        assertEquals(35, readState.refPosition());
        assertEquals(34, readState.readIndex());
        assertEquals(M, readState.operator());

        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(31, readState.refPosition());
        assertEquals(29, readState.readIndex());
        assertEquals(I, readState.operator());

        // move through the insert
        readState.moveNext();

        assertEquals(31, readState.refPosition());
        assertEquals(28, readState.readIndex());
        assertEquals(I, readState.operator());

        readState.moveNext();
        readState.moveNext();
        readState.moveNext();
        readState.moveNext();

        assertEquals(30, readState.refPosition());
        assertEquals(24, readState.readIndex());
        assertEquals(M, readState.operator());

        // move to start of delete
        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(25, readState.refPosition());
        assertEquals(20, readState.readIndex());
        assertEquals(D, readState.operator());

        // move past delete
        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(20, readState.refPosition());
        assertEquals(19, readState.readIndex());
        assertEquals(M, readState.operator());

        // move to last aligned base
        for(int i = 0; i < 9; ++i)
        {
            readState.moveNext();
        }

        assertEquals(11, readState.refPosition());
        assertEquals(10, readState.readIndex());
        assertEquals(M, readState.operator());

        // stop at final soft-clipping
        readState.moveNext();
        assertTrue(readState.exhausted());
    }

    @Test
    public void testHardClippedReads()
    {
        // index 0-14, ref pos 1-15
        // index DEL, ref pos 16-21
        // index 15-46, ref pos 22-52
        String cigar = "105H15M6D31M";
        String readBases = REF_BASES_200.substring(1, 16) + REF_BASES_200.substring(22, 53);

        Read read = createRead(READ_ID_GENERATOR.nextId(), 1, readBases, cigar);
        int readJunctionIndex = 15;
        boolean isForwardJunction = false;

        RefReadParseState readState = new RefReadParseState(
                isForwardJunction, read, readJunctionIndex, readRefBaseLength(read, readJunctionIndex, isForwardJunction));

        assertEquals(22, readState.refPosition());
        assertEquals(15, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(22), (char)readState.currentBase());

        readState.moveNext();

        assertEquals(23, readState.refPosition());
        assertEquals(16, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(23), (char)readState.currentBase());

        // index 0-30, ref pos 1-31
        // index DEL, ref pos 32-37
        // index 31-45, ref pos 38-52
        cigar = "31M6D15M105H";
        readBases = REF_BASES_200.substring(1, 32) + REF_BASES_200.substring(38, 53);

        read = createRead(READ_ID_GENERATOR.nextId(), 1, readBases, cigar);
        readJunctionIndex = 31;
        isForwardJunction = true;

        readState = new RefReadParseState(
                isForwardJunction, read, readJunctionIndex, readRefBaseLength(read, readJunctionIndex, isForwardJunction));

        assertEquals(38, readState.refPosition());
        assertEquals(31, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(38), (char)readState.currentBase());

        readState.moveNext();

        assertEquals(37, readState.refPosition());
        assertEquals(31, readState.readIndex());
        assertEquals(D, readState.operator());

        for(int i = 0; i < 6; ++i)
        {
            readState.moveNext();
        }

        assertEquals(31, readState.refPosition());
        assertEquals(30, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(31), (char)readState.currentBase());
    }

    @Test
    public void testForwardRefBaseSequences()
    {
        String extBases = REF_BASES_200.substring(100, 140);
        byte[] extBaseQuals = buildDefaultBaseQuals(extBases.length());

        Junction junction = new Junction(CHR_1, 100, FORWARD);

        // first a basic assembly with all reads agreeing
        String readBases = REF_BASES_200.substring(81, 101) + extBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 81, readBases, makeCigarString(readBases, 0, extBases.length()));

        readBases = REF_BASES_200.substring(51, 101) + extBases; // extends but stops at soft-clip
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 71, readBases, makeCigarString(readBases, 20, extBases.length()));

        readBases = REF_BASES_200.substring(61, 101) + extBases;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 61, readBases, makeCigarString(readBases, 0, extBases.length()));

        List<Read> reads = List.of(read1, read2, read3);

        List<SupportRead> supportReads = reads.stream()
                .map(x -> new SupportRead(x, SupportType.JUNCTION, 0, 0, 0)).collect(Collectors.toList());

        JunctionAssembly assembly = new JunctionAssembly(junction, extBases.getBytes(), extBaseQuals, supportReads, Collections.emptyList());

        RefBaseSeqBuilder refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);

        String refBases = REF_BASES_200.substring(61, 101);
        assertEquals(refBases, refBaseSeqBuilder.refBaseSequence());
        assertEquals(61, refBaseSeqBuilder.refBasePosition());
        assertEquals(40, refBaseSeqBuilder.refBaseLength());
        assertEquals("40M", refBaseSeqBuilder.cigarStr());

        assertTrue(refBaseSeqBuilder.reads().stream().allMatch(x -> x.mismatches() == 0));

        // now with reads with disagreeing bases
        String readRefBases = "AAACCCGGGTTTACGTAACCGGTTACGTAA";
        //                     012345678901234567890123456789

        readBases = readRefBases + extBases;
        read1 = createRead(READ_ID_GENERATOR.nextId(), 71, readBases, makeCigarString(readBases, 0, extBases.length()));

        readBases = readRefBases.substring(0, 5) + "A" + readRefBases.substring(6, 16) + "T" + readRefBases.substring(17) + extBases;
        read2 = createRead(READ_ID_GENERATOR.nextId(), 71, readBases, makeCigarString(readBases, 0, extBases.length()));

        readBases = readRefBases.substring(0, 5) + "A" + readRefBases.substring(6, 26) + "T" + readRefBases.substring(27) + extBases;
        read3 = createRead(READ_ID_GENERATOR.nextId(), 71, readBases, makeCigarString(readBases, 0, extBases.length()));

        reads = List.of(read1, read2, read3);

        supportReads = reads.stream()
                .map(x -> new SupportRead(x, SupportType.JUNCTION, 0, 0, 0)).collect(Collectors.toList());

        String juncBases = readRefBases.substring(29) + extBases;
        assembly = new JunctionAssembly(junction, juncBases.getBytes(), extBaseQuals, supportReads, Collections.emptyList());

        refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);

        String refSeqBases = readRefBases.substring(0, 5) + "A" + readRefBases.substring(6);
        assertEquals(refSeqBases, refBaseSeqBuilder.refBaseSequence());
        assertEquals(71, refBaseSeqBuilder.refBasePosition());
        assertEquals(30, refBaseSeqBuilder.refBaseLength());
        assertEquals("30M", refBaseSeqBuilder.cigarStr());

        assertEquals(3, refBaseSeqBuilder.reads().stream().filter(x -> x.mismatches() == 1).count());

        // with a mix of INDELs

        //            71       80          90               100
        //            1234567890 12345 67890       1234567890
        // refBases1: AAACCCGGGT       TAACC TTTTT GGTTACGTAA"
        // refBases2: AAACCCGGGT TTACG TAACC TTTTT GGTTACGTAA"
        // refBases3: AAACCCGGGT       TAACC       GGTTACGTAA"
        String insert = "TT";
        String refBases1 = readRefBases.substring(0, 10) + readRefBases.substring(15, 20) + insert + readRefBases.substring(20, 30);
        read1 = createRead(READ_ID_GENERATOR.nextId(), 71, refBases1 + extBases, "10M5D5M2I10M40S");

        // has the insert but not the delete
        String refBases2 = readRefBases.substring(0, 20) + insert + readRefBases.substring(20, 30);
        read2 = createRead(READ_ID_GENERATOR.nextId(), 71, refBases2 + extBases, "20M2I10M40S");

        // has the delete but not the insert
        String refBases3 = readRefBases.substring(0, 10) + readRefBases.substring(15, 30);
        read3 = createRead(READ_ID_GENERATOR.nextId(), 71, refBases3 + extBases, "10M5D15M40S");

        reads = List.of(read1, read2, read3);

        supportReads = reads.stream()
                .map(x -> new SupportRead(x, SupportType.JUNCTION, 0, 0, 0)).collect(Collectors.toList());

        juncBases = readRefBases.substring(29) + extBases;
        assembly = new JunctionAssembly(junction, juncBases.getBytes(), extBaseQuals, supportReads, Collections.emptyList());

        refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);

        refSeqBases = readRefBases.substring(0, 10) + readRefBases.substring(15, 20) + insert + readRefBases.substring(20, 30);
        assertEquals(refSeqBases, refBaseSeqBuilder.refBaseSequence());
        assertEquals(71, refBaseSeqBuilder.refBasePosition());
        assertEquals(27, refBaseSeqBuilder.refBaseLength());
        assertEquals("10M5D5M2I10M", refBaseSeqBuilder.cigarStr());

        assertEquals(0, getReadMismatchCount(refBaseSeqBuilder, read1));
        assertEquals(5, getReadIndelMismatchCount(refBaseSeqBuilder, read2));
        assertEquals(2, getReadIndelMismatchCount(refBaseSeqBuilder, read3));
    }

    @Test
    public void testReverseRefBaseSequences()
    {
        // matches the test above
        String extBases = REF_BASES_200.substring(20, 40);

        Junction junction = new Junction(CHR_1, 100, REVERSE);

        String readBases = extBases + REF_BASES_200.substring(100, 120);
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, makeCigarString(readBases, extBases.length(), 0));

        readBases = extBases + REF_BASES_200.substring(100, 150);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, makeCigarString(readBases, extBases.length(), 20));

        readBases = extBases + REF_BASES_200.substring(100, 140);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, makeCigarString(readBases, extBases.length(), 0));

        List<Read> reads = List.of(read1, read2, read3);

        List<SupportRead> supportReads = reads.stream()
                .map(x -> new SupportRead(x, SupportType.JUNCTION, 0, 0, 0)).collect(Collectors.toList());

        String juncBases = extBases + REF_BASES_200.substring(100, 101);
        byte[] juncBaseQuals = buildDefaultBaseQuals(juncBases.length());
        JunctionAssembly assembly = new JunctionAssembly(junction, juncBases.getBytes(), juncBaseQuals, supportReads, Collections.emptyList());

        RefBaseSeqBuilder refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);

        String refBases = REF_BASES_200.substring(100, 140);
        assertEquals(refBases, refBaseSeqBuilder.refBaseSequence());
        assertEquals(139, refBaseSeqBuilder.refBasePosition());
        assertEquals(40, refBaseSeqBuilder.refBaseLength());
        assertEquals("40M", refBaseSeqBuilder.cigarStr());

        assertTrue(refBaseSeqBuilder.reads().stream().allMatch(x -> x.mismatches() == 0));

        String readRefBases = "AAACCCGGGTTTACGTAACCGGTTACGTAA";
        //                     012345678901234567890123456789

        // with a mix of INDELs

        //            100        110               120       130
        //            0123456789 01234 56789       01234567890
        // refBases1: AAACCCGGGT       TAACC TTTTT GGTTACGTAA
        // refBases2: AAACCCGGGT TTACG TAACC TTTTT GGTTACGTAA
        // refBases3: AAACCCGGGT       TAACC       GGTTACGTAA
        String insert = "TTTTT";
        String refBases1 = readRefBases.substring(0, 10) + readRefBases.substring(15, 20) + insert + readRefBases.substring(20, 30);
        read1 = createRead(READ_ID_GENERATOR.nextId(), 100, extBases + refBases1, "20S10M5D5M5I10M");

        // has the insert but not the delete
        String refBases2 = readRefBases.substring(0, 20) + insert + readRefBases.substring(20, 30);
        read2 = createRead(READ_ID_GENERATOR.nextId(), 100, extBases + refBases2, "20S20M5I10M");

        // has the delete but not the insert
        String refBases3 = readRefBases.substring(0, 10) + readRefBases.substring(15, 30);
        read3 = createRead(READ_ID_GENERATOR.nextId(), 100, extBases + refBases3, "20S10M5D15M");

        reads = List.of(read1, read2, read3);

        supportReads = reads.stream()
                .map(x -> new SupportRead(x, SupportType.JUNCTION, 0, 0, 0)).collect(Collectors.toList());

        juncBases = extBases + readRefBases.substring(0, 1);
        juncBaseQuals = buildDefaultBaseQuals(juncBases.length());

        assembly = new JunctionAssembly(junction, juncBases.getBytes(), juncBaseQuals, supportReads, Collections.emptyList());

        refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);

        String refSeqBases = readRefBases.substring(0, 10) + readRefBases.substring(15, 20) + insert + readRefBases.substring(20, 30);
        assertEquals(refSeqBases, refBaseSeqBuilder.refBaseSequence());
        assertEquals(129, refBaseSeqBuilder.refBasePosition());
        assertEquals(30, refBaseSeqBuilder.refBaseLength());
        assertEquals("10M5D5M5I10M", refBaseSeqBuilder.cigarStr());

        assertEquals(0, getReadMismatchCount(refBaseSeqBuilder, read1));
        assertEquals(5, getReadIndelMismatchCount(refBaseSeqBuilder, read2));
        assertEquals(3, getReadIndelMismatchCount(refBaseSeqBuilder, read3));
    }

    @Test
    public void testSupportingIndelTolerances()
    {
        String extBases = REF_BASES_200.substring(100, 140);
        byte[] extBaseQuals = buildDefaultBaseQuals(extBases.length());

        Junction junction = new Junction(CHR_1, 100, FORWARD);

        // the 2 reads cover the same ref bases (70-100) but have different cigar representations
        String readBases = REF_BASES_200.substring(71, 101) + extBases;

        List<SupportRead> supportReads = Lists.newArrayList();

        String cigar = makeCigarString(readBases, 0, extBases.length());
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 71, readBases, cigar);
        Read read2 = cloneRead(read1, READ_ID_GENERATOR.nextId());

        supportReads.add(new SupportRead(read1, SupportType.JUNCTION, 30, 0, 0));
        supportReads.add(new SupportRead(read2, SupportType.JUNCTION, 30, 0, 0));

        String readBases2 = REF_BASES_200.substring(61, 101) + extBases.substring(0, 30);

        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 61, readBases2, "30M10I30M");
        Read read4 = cloneRead(read3, READ_ID_GENERATOR.nextId());
        Read read5 = cloneRead(read3, READ_ID_GENERATOR.nextId());

        calcIndelInferredUnclippedPositions(read3);
        calcIndelInferredUnclippedPositions(read4);
        calcIndelInferredUnclippedPositions(read5);

        supportReads.add(new SupportRead(read3, SupportType.JUNCTION, 30, 0, 0));
        supportReads.add(new SupportRead(read4, SupportType.JUNCTION, 30, 0, 0));
        supportReads.add(new SupportRead(read5, SupportType.JUNCTION, 30, 0, 0));

        JunctionAssembly assembly = new JunctionAssembly(junction, extBases.getBytes(), extBaseQuals, supportReads, Collections.emptyList());

        RefBaseSeqBuilder refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);

        String refBases = REF_BASES_200.substring(61, 101);

        assertEquals(refBases, refBaseSeqBuilder.refBaseSequence());
        assertEquals(61, refBaseSeqBuilder.refBasePosition());
        assertEquals(40, refBaseSeqBuilder.refBaseLength());
        assertEquals("30M9I1M", refBaseSeqBuilder.cigarStr());

        assertTrue(refBaseSeqBuilder.reads().stream().allMatch(x -> x.mismatches() == 0));
    }

    private static int getReadMismatchCount(final RefBaseSeqBuilder refBaseSeqBuilder, final Read read)
    {
        RefReadParseState readState = refBaseSeqBuilder.reads().stream().filter(x -> x.read() == read).findFirst().orElse(null);
        return readState.mismatches();
    }

    private static int getReadIndelMismatchCount(final RefBaseSeqBuilder refBaseSeqBuilder, final Read read)
    {
        RefReadParseState readState = refBaseSeqBuilder.reads().stream().filter(x -> x.read() == read).findFirst().orElse(null);
        return readState.indelMismatches();
    }
}
