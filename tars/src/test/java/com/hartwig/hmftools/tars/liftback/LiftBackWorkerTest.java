package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.LiftBackWorker.sanitizeForOutput;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.mappedRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.pairedUnmappedRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.repeatedBase;

import static org.junit.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

// Write-boundary guard: htsjdk and REDUX both reject a zero-length CIGAR element or a SEQ length disagreeing with the CIGAR.
public class LiftBackWorkerTest
{
    private static SAMRecord mappedRecordWithSeq(final Cigar cigar, final int seqLength)
    {
        SAMRecord record = mappedRecord("read1", CHR_1, 100, cigar.toString());
        record.setReadBases(repeatedBase(seqLength, 'A'));
        return record;
    }

    @Test
    public void testSanitizeStripsZeroLengthCigarElement()
    {
        Cigar withZero = new Cigar(List.of(new CigarElement(0, CigarOperator.M), new CigarElement(50, CigarOperator.M)));
        SAMRecord record = mappedRecordWithSeq(withZero, 50);

        sanitizeForOutput(record);

        assertEquals("50M", record.getCigarString());
    }

    @Test
    public void testSanitizeSeqCigarMismatchReplacedWithPlaceholder()
    {
        // a failed-lift supplementary mirrored onto its primary's 100M cigar while carrying only 50 hard-clipped bases
        SAMRecord record = mappedRecordWithSeq(new Cigar(List.of(new CigarElement(100, CigarOperator.M))), 50);

        sanitizeForOutput(record);

        assertEquals("50M", record.getCigarString());
        assertEquals(50, record.getCigar().getReadLength());
    }

    @Test
    public void testSanitizeValidRecordUnchanged()
    {
        SAMRecord record = mappedRecordWithSeq(new Cigar(List.of(
                new CigarElement(10, CigarOperator.S), new CigarElement(40, CigarOperator.M))), 50);

        sanitizeForOutput(record);

        assertEquals("10S40M", record.getCigarString());
    }

    @Test
    public void testSanitizeUnmappedNoOp()
    {
        SAMRecord record = pairedUnmappedRecord("read1", true);
        record.setCigarString("50M");

        sanitizeForOutput(record);

        assertEquals("50M", record.getCigarString());
    }
}
