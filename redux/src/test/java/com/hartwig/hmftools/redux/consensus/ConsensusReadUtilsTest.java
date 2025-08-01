package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.redux.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.redux.TestUtils.createFragmentCoords;
import static com.hartwig.hmftools.redux.consensus.TemplateReads.selectTemplateRead;
import static com.hartwig.hmftools.redux.umi.UmiConfig.READ_ID_DELIM_STR;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;
import com.hartwig.hmftools.redux.consensus.ReadParseState;

import org.junit.jupiter.api.Test;

import htsjdk.samtools.SAMRecord;

public class ConsensusReadUtilsTest
{
    @Test
    public void testConsensusReadId()
    {
        String readIdFixed = "ABAB:8:SAMPLE:2:222:12345";

        SAMRecord read1 = createSamRecord(
                readIdFixed + READ_ID_DELIM_STR + "READ_01", 100, TEST_READ_BASES, TEST_READ_CIGAR, TEST_READ_CIGAR);

        String consensusReadId = ConsensusReads.formConsensusReadId(read1, null);

        assertEquals("ABAB:8:SAMPLE:2:222:12345:CNS_READ_01", consensusReadId);

        String unmiId = "ACGTG_ATTGC";

        read1 = createSamRecord(
                readIdFixed + READ_ID_DELIM_STR + "READ_01" + READ_ID_DELIM_STR + unmiId,
                100, TEST_READ_BASES, TEST_READ_CIGAR, TEST_READ_CIGAR);

        consensusReadId = ConsensusReads.formConsensusReadId(read1, unmiId);

        assertEquals("ABAB:8:SAMPLE:2:222:12345:READ_01:CNS_" + unmiId, consensusReadId);
    }

    @Test
    public void testTemplateReadSelection()
    {
        int posStart = 11;
        String consensusBases = REF_BASES.substring(posStart, 21);

        // all agreeing
        SAMRecord read1 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, TEST_READ_CIGAR);
        SAMRecord read2 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, TEST_READ_CIGAR);
        SAMRecord read3 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, TEST_READ_CIGAR);

        FragmentCoords fragmentCoords = createFragmentCoords(read1);

        SAMRecord templateRead = selectTemplateRead(List.of(read3, read2, read1), fragmentCoords);
        assertEquals(read1, templateRead);

        // select on lower cigar
        read1 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, "10S90M", TEST_READ_CIGAR);
        read2 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, TEST_READ_CIGAR);
        read3 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, TEST_READ_CIGAR);

        fragmentCoords = createFragmentCoords(read3);
        templateRead = selectTemplateRead(List.of(read3, read2, read1), fragmentCoords);
        assertEquals(read2, templateRead);

        // select on mate cigar
        read1 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, "90M10S");
        read2 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, "10S90M");
        read3 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, TEST_READ_CIGAR);
        SAMRecord read4 = createSamRecord(READ_ID_GEN.nextId(), posStart, consensusBases, TEST_READ_CIGAR, TEST_READ_CIGAR);

        fragmentCoords = createFragmentCoords(read4);
        templateRead = selectTemplateRead(List.of(read4, read3, read2, read1), fragmentCoords);
        assertEquals(read3, templateRead);
    }

    private static SAMRecord createSamRecord(
            final String readId, int readStart, final String readBases, final String cigar, final String mateCigar)
    {
        return SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, 5000, false, false, null,
                true, mateCigar);
    }

    @Test
    public void testReadParseState()
    {
        String bases = "AGGCGGA";
        String indelCigar = "1S2M1I2M1S";

        SAMRecord read1 = createSamRecord(TEST_READ_ID, 100, bases, indelCigar, TEST_READ_CIGAR);

        ReadParseState readState = new ReadParseState(read1, true);
        assertEquals((byte)'A', readState.currentBase());
        assertEquals(DEFAULT_QUAL, readState.currentBaseQual());
        assertEquals(S, readState.elementType());
        assertEquals(1, readState.elementLength());

        readState.moveNextBase();
        assertEquals((byte)'G', readState.currentBase());
        assertEquals(M, readState.elementType());
        assertEquals(2, readState.elementLength());

        readState.moveNextBase();
        readState.moveNextBase();
        assertEquals((byte)'C', readState.currentBase());
        assertEquals(I, readState.elementType());
        assertEquals(1, readState.elementLength());

        readState.moveNextBase();
        readState.moveNextBase();
        assertFalse(readState.exhausted());
        assertEquals((byte)'G', readState.currentBase());
        assertEquals(M, readState.elementType());

        readState.moveNextBase();
        assertFalse(readState.exhausted());
        assertEquals((byte)'A', readState.currentBase());
        assertEquals(S, readState.elementType());

        readState.moveNextBase();
        assertTrue(readState.exhausted());

        // and in reverse
        readState = new ReadParseState(read1, false);
        assertEquals((byte)'A', readState.currentBase());
        assertEquals(DEFAULT_QUAL, readState.currentBaseQual());
        assertEquals(S, readState.elementType());
        assertEquals(1, readState.elementLength());

        readState.moveNextBase();
        readState.moveNextBase();
        readState.moveNextBase();

        assertEquals((byte)'C', readState.currentBase());
        assertEquals(I, readState.elementType());
        assertEquals(1, readState.elementLength());

        readState.moveNextBase();
        readState.moveNextBase();
        readState.moveNextBase();
        assertFalse(readState.exhausted());
        assertEquals((byte)'A', readState.currentBase());
        assertEquals(S, readState.elementType());

        readState.moveNextBase();
        assertTrue(readState.exhausted());

        // with a delete
        bases = "ACGT";
        indelCigar = "2M3D2M";

        read1 = createSamRecord(TEST_READ_ID, 100, bases, indelCigar, TEST_READ_CIGAR);

        readState = new ReadParseState(read1, true);

        assertEquals((byte)'A', readState.currentBase());
        assertEquals(M, readState.elementType());

        readState.moveNextBase();
        assertEquals((byte)'C', readState.currentBase());
        assertEquals(M, readState.elementType());

        readState.moveNextBase();
        assertEquals((byte)'C', readState.currentBase());
        assertEquals(D, readState.elementType());

        readState.moveNextBase();
        assertEquals((byte)'C', readState.currentBase());
        assertEquals(D, readState.elementType());

        readState.moveNextBase();
        assertEquals((byte)'C', readState.currentBase());
        assertEquals(D, readState.elementType());

        readState.moveNextBase();
        assertEquals((byte)'G', readState.currentBase());
        assertEquals(M, readState.elementType());

        readState.moveNextBase();
        assertEquals((byte)'T', readState.currentBase());
        assertEquals(M, readState.elementType());

        readState.moveNextBase();
        assertTrue(readState.exhausted());
    }
}
