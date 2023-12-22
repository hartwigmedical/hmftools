package com.hartwig.hmftools.svassembly.sam;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;

import com.hartwig.hmftools.svassembly.models.Record;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class PolyGTrimmerTest
{
    private static Record record(final String bases, final boolean isPositiveStrand)
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

        return new Record(record);
    }

    @Test
    public void trimsGsForward()
    {
        final PolyGTrimmer trimmer = new PolyGTrimmer(3);

        final var input = record("CCCCCAAAAAGGGGG", true);
        final var trimmed = trimmer.trimPolyG(input);

        assertThat(trimmed).isNotSameAs(input);
        assertThat(trimmed.getBasesString()).isEqualTo("CCCCCAAAAA");
        assertThat(trimmed.getAlignmentStart()).isEqualTo(100);
        assertThat(trimmed.getAlignmentEnd()).isEqualTo(109);
        assertThat(trimmed.isPositiveStrand()).isEqualTo(true);
        assertThat(trimmed.getCigar().toString()).isEqualTo("10M");
    }

    @Test
    public void trimsCsBackward()
    {
        final PolyGTrimmer trimmer = new PolyGTrimmer(3);

        final var input = record("CCCCCAAAAAGGGGG", false);
        final var trimmed = trimmer.trimPolyG(input);

        assertThat(trimmed).isNotSameAs(input);
        assertThat(trimmed.getBasesString()).isEqualTo("AAAAAGGGGG");
        assertThat(trimmed.getAlignmentStart()).isEqualTo(105);
        assertThat(trimmed.getAlignmentEnd()).isEqualTo(114);
        assertThat(trimmed.isPositiveStrand()).isEqualTo(false);
        assertThat(trimmed.getCigar().toString()).isEqualTo("10M");
    }

    @Test
    public void dontTrimGsBelowThreshold()
    {
        final PolyGTrimmer trimmer = new PolyGTrimmer(10);

        final var input = record("CCCCCAAAAAGGGGG", true);
        final var trimmed = trimmer.trimPolyG(input);

        assertThat(trimmed).isSameAs(input);
        assertThat(trimmed.getBasesString()).isEqualTo("CCCCCAAAAAGGGGG"); // Check we didn't change anything
    }
}