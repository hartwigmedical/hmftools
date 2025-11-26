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
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.BASE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

import org.junit.Test;

public class RefBaseSequenceTest
{
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

        assertTrue(refBaseSeqBuilder.reads().stream().allMatch(x -> x.mismatchCount(true) == 0));

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

        assertEquals(3, refBaseSeqBuilder.reads().stream().filter(x -> x.mismatchCount(true) == 1).count());

        // with a mix of INDELs, all marked as low-qual

        //            71       80          90            100
        //            1234567890 12345 67890    1234567890
        // refBases1: AAACCCGGGT       TAACC TT GGTTACGTAA"
        // refBases2: AAACCCGGGT TTACG TAACC TT GGTTACGTAA"
        // refBases3: AAACCCGGGT       TAACC    GGTTACGTAA"
        String insert = "TT";
        String refBases1 = readRefBases.substring(0, 10) + readRefBases.substring(15, 20) + insert + readRefBases.substring(20, 30);
        read1 = createRead(READ_ID_GENERATOR.nextId(), 71, refBases1 + extBases, "10M5D5M2I10M40S");

        // has the insert but not the delete
        String refBases2 = readRefBases.substring(0, 20) + insert + readRefBases.substring(20, 30);
        read2 = createRead(READ_ID_GENERATOR.nextId(), 71, refBases2 + extBases, "20M2I10M40S");
        setBaseQuals(read2, 10, 14, LOW_QUAL_BASE);

        // has the delete but not the insert
        String refBases3 = readRefBases.substring(0, 10) + readRefBases.substring(15, 30);
        read3 = createRead(READ_ID_GENERATOR.nextId(), 71, refBases3 + extBases, "10M5D15M40S");
        setBaseQuals(read3, 14, 15, LOW_QUAL_BASE);

        reads = List.of(read1, read2, read3);

        supportReads = reads.stream()
                .map(x -> new SupportRead(x, SupportType.JUNCTION, 0, 0, 0)).collect(Collectors.toList());

        juncBases = readRefBases.substring(29) + extBases;
        assembly = new JunctionAssembly(junction, juncBases.getBytes(), extBaseQuals, supportReads, Collections.emptyList());

        refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);

        refSeqBases = readRefBases.substring(0, 10) + readRefBases.substring(15, 20) + insert + readRefBases.substring(20, 30);
        assertEquals("10M5D5M2I10M", refBaseSeqBuilder.cigarStr());
        assertEquals(refSeqBases, refBaseSeqBuilder.refBaseSequence());
        assertEquals(71, refBaseSeqBuilder.refBasePosition());
        assertEquals(27, refBaseSeqBuilder.refBaseLength());

        assertEquals(0, getReadMismatchCount(refBaseSeqBuilder, read1));
        assertEquals(5, getReadMismatchCount(refBaseSeqBuilder, read2));
        assertEquals(2, getReadMismatchCount(refBaseSeqBuilder, read3));
    }

    private final byte LOW_QUAL_BASE = BaseQualAdjustment.LOW_BASE_QUAL_THRESHOLD - 1; // until isLowQual is corrected to include 26

    private static void setBaseQuals(final Read read, int readStart, int readEnd, byte qual)
    {
        byte[] readQuals = read.getBaseQuality();
        for(int i = readStart; i <= readEnd; ++i)
        {
            readQuals[i] = qual;
        }
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

        assertTrue(refBaseSeqBuilder.reads().stream().allMatch(x -> x.mismatchCount(true) == 0));

        String readRefBases = "AAACCCGGGTTTACGTAACCGGTTACGTAA";
        //                     012345678901234567890123456789

        // with a mix of INDELs

        //            100        110            120       130
        //            0123456789 01234 56789    01234567890
        // refBases1: AAACCCGGGT       TAACC TT GGTTACGTAA
        // refBases2: AAACCCGGGT TTACG TAACC TT GGTTACGTAA
        // refBases3: AAACCCGGGT       TAACC    GGTTACGTAA
        String insert = "TT";
        String refBases1 = readRefBases.substring(0, 10) + readRefBases.substring(15, 20) + insert + readRefBases.substring(20, 30);
        read1 = createRead(READ_ID_GENERATOR.nextId(), 100, extBases + refBases1, "20S10M5D5M2I10M");

        // has the insert but not the delete
        String refBases2 = readRefBases.substring(0, 20) + insert + readRefBases.substring(20, 30);
        read2 = createRead(READ_ID_GENERATOR.nextId(), 100, extBases + refBases2, "20S20M2I10M");
        setBaseQuals(read2, 30, 34, LOW_QUAL_BASE);

        // has the delete but not the insert
        String refBases3 = readRefBases.substring(0, 10) + readRefBases.substring(15, 30);
        read3 = createRead(READ_ID_GENERATOR.nextId(), 100, extBases + refBases3, "20S10M5D15M");
        setBaseQuals(read3, 34, 35, LOW_QUAL_BASE);

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
        assertEquals(27, refBaseSeqBuilder.refBaseLength());
        assertEquals("10M5D5M2I10M", refBaseSeqBuilder.cigarStr());

        assertEquals(0, getReadMismatchCount(refBaseSeqBuilder, read1));
        assertEquals(5, getReadMismatchCount(refBaseSeqBuilder, read2));
        assertEquals(2, getReadMismatchCount(refBaseSeqBuilder, read3));
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
        setBaseQuals(read1, 20, 29, LOW_QUAL_BASE);
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

        assertTrue(refBaseSeqBuilder.reads().stream().allMatch(x -> x.mismatchCount(true) == 0));
    }

    private static int getReadMismatchCount(final RefBaseSeqBuilder refBaseSeqBuilder, final Read read)
    {
        ReadParseState readState = refBaseSeqBuilder.reads().stream().filter(x -> x.read() == read).findFirst().orElse(null);
        return (int)readState.mismatches().stream().filter(x -> x.Type == BASE).count();
    }

    private static int getReadIndelMismatchCount(final RefBaseSeqBuilder refBaseSeqBuilder, final Read read)
    {
        ReadParseState readState = refBaseSeqBuilder.reads().stream().filter(x -> x.read() == read).findFirst().orElse(null);
        return (int)readState.mismatches().stream().filter(x -> x.Type != BASE).count();
    }
}
