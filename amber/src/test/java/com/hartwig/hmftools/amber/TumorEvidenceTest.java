package com.hartwig.hmftools.amber;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.bam.SamRecordUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class TumorEvidenceTest
{
    @Test
    public void useQualityOfBaseAfterDel()
    {
        int minQuality = SamRecordUtils.getBaseQuality('J');

        SAMRecord lowBaseQualDel = buildSamRecord(1000, "1M1D1M", "CT", "FI");
        lowBaseQualDel.setMappingQuality(60);

        SAMRecord lowMapQualRead = buildSamRecord(1000, "1M1D1M", "CT", "FJ");
        lowMapQualRead.setMappingQuality(40);

        final PositionEvidence baseDepth = new PositionEvidence("5", 1001, "A", "T");

        PositionEvidenceChecker evidenceChecker = new PositionEvidenceChecker(50, minQuality);

        evidenceChecker.addEvidence(baseDepth, lowBaseQualDel);
        assertEquals(1, baseDepth.ReadDepth);
        assertEquals(1, baseDepth.BaseQualFiltered);

        evidenceChecker.addEvidence(baseDepth, lowMapQualRead);
        assertEquals(2, baseDepth.ReadDepth);
        assertEquals(1, baseDepth.MapQualFiltered);
    }

    private SAMRecord buildSamRecord(
            final int alignmentStart, final String cigar, final String readString, final String qualities)
    {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        return record;
    }
}
