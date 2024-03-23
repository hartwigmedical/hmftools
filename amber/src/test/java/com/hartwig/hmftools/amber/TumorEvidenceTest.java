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

        final SAMRecord lowQualDel = buildSamRecord(1000, "1M1D1M", "CT", "FI");
        final SAMRecord highQualDel = buildSamRecord(1000, "1M1D1M", "CT", "FJ");

        final PositionEvidence baseDepth = new PositionEvidence("5", 1001, "A", "T");

        PositionEvidenceChecker evidenceChecker = new PositionEvidenceChecker(minQuality);

        evidenceChecker.addEvidence(baseDepth, lowQualDel);
        assertEquals(0, baseDepth.ReadDepth);

        evidenceChecker.addEvidence(baseDepth, highQualDel);
        assertEquals(1, baseDepth.ReadDepth);
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
