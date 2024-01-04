package com.hartwig.hmftools.esvee.read;

import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class PolyGTrimmerTest
{
    private static Read record(final String bases, final boolean isPositiveStrand)
    {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("1", 10_000));

        final SAMRecord record = new SAMRecord(header);
        record.setReadBases(bases.getBytes());
        record.setBaseQualities(new byte[bases.length()]);
        record.setReadPairedFlag(true);
        record.setReadNegativeStrandFlag(!isPositiveStrand);
        record.setReferenceName("1");
        record.setAlignmentStart(100);
        record.setCigar(new Cigar(List.of(new CigarElement(bases.length(), CigarOperator.M))));

        return new Read(record);
    }

    /* CHASHA FIXME
    @Test
    public void trimsGsForward()
    {
        final PolyGTrimmer trimmer = new PolyGTrimmer(3);

        final var input = record("CCCCCAAAAAGGGGG", true);
        final var trimmed = trimmer.trimPolyG(input);

        assertTrue(trimmed).isNotSameAs(input);
        assertTrue(trimmed.getBasesString()).isEqualTo("CCCCCAAAAA");
        assertTrue(trimmed.getAlignmentStart()).isEqualTo(100);
        assertTrue(trimmed.getAlignmentEnd()).isEqualTo(109);
        assertTrue(trimmed.isPositiveStrand()).isEqualTo(true);
        assertTrue(trimmed.getCigar().toString()).isEqualTo("10M");
    }

    @Test
    public void trimsCsBackward()
    {
        final PolyGTrimmer trimmer = new PolyGTrimmer(3);

        final var input = record("CCCCCAAAAAGGGGG", false);
        final var trimmed = trimmer.trimPolyG(input);

        assertTrue(trimmed).isNotSameAs(input);
        assertTrue(trimmed.getBasesString()).isEqualTo("AAAAAGGGGG");
        assertTrue(trimmed.getAlignmentStart()).isEqualTo(105);
        assertTrue(trimmed.getAlignmentEnd()).isEqualTo(114);
        assertTrue(trimmed.isPositiveStrand()).isEqualTo(false);
        assertTrue(trimmed.getCigar().toString()).isEqualTo("10M");
    }

    @Test
    public void dontTrimGsBelowThreshold()
    {
        final PolyGTrimmer trimmer = new PolyGTrimmer(10);

        final var input = record("CCCCCAAAAAGGGGG", true);
        final var trimmed = trimmer.trimPolyG(input);

        assertTrue(trimmed).isSameAs(input);
        assertTrue(trimmed.getBasesString()).isEqualTo("CCCCCAAAAAGGGGG"); // Check we didn't change anything
    }
    */
}