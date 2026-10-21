package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MIN_BASE_QUALITY;
import static com.hartwig.hmftools.amber.AmberConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.bam.SamRecordUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class TumorEvidenceTest
{
    private static final String READ_ID = "READ_01";

    @Test
    public void useQualityOfBaseAfterDel()
    {
        int minQuality = SamRecordUtils.getBaseQuality('J');

        SAMRecord lowBaseQualDel = createSamRecord(
                READ_ID, CHR_1, 1000, "CT", "1M1D1M", CHR_1, 1200, false,
                false, null);
        lowBaseQualDel.setBaseQualityString("FI");

        SAMRecord lowMapQualRead = createSamRecord(
                READ_ID, CHR_1, 1000, "CT", "1M1D1M", CHR_1, 1200, false,
                false, null);
        lowMapQualRead.setMappingQuality(40);

        PositionEvidence baseDepth = new PositionEvidence("5", 1001, "A", "T");

        PositionEvidenceChecker evidenceChecker = new PositionEvidenceChecker(50, minQuality);

        ReadCigarState state = new ReadCigarState(
                lowBaseQualDel.getCigar().getCigarElements(), lowBaseQualDel.getAlignmentStart(), lowBaseQualDel.getReadLength());
        state.moveToPosition(1001);

        evidenceChecker.addEvidence(baseDepth, lowBaseQualDel, state.ReadIndex, state.isIndel());
        assertEquals(1, baseDepth.ReadDepth);
        assertEquals(1, baseDepth.BaseQualFiltered);

        state = new ReadCigarState(
                lowMapQualRead.getCigar().getCigarElements(), lowMapQualRead.getAlignmentStart(), lowMapQualRead.getReadLength());
        state.moveToPosition(1001);

        evidenceChecker.addEvidence(baseDepth, lowMapQualRead, state.ReadIndex, state.isIndel());
        assertEquals(2, baseDepth.ReadDepth);
        assertEquals(1, baseDepth.MapQualFiltered);
    }

    @Test
    public void testFragmentSync()
    {
        PositionEvidence posEvidence1 = new PositionEvidence(CHR_1, 1020, "A", "C");
        PositionEvidence posEvidence2 = new PositionEvidence(CHR_1, 1025, "A", "C");

        RegionTask regionTask = new RegionTask(
                new PositionEvidenceChecker(DEFAULT_MIN_MAPPING_QUALITY, DEFAULT_MIN_BASE_QUALITY), CHR_1, posEvidence1);
        regionTask.addPosition(posEvidence2);

        String readBases = "A".repeat(30);

        // both reads overlap the same 2 positions but are only counted once
        SAMRecord read1 = createSamRecord(
                READ_ID, CHR_1, 1000, readBases, "30M", CHR_1, 1010, false,
                false, null, true, "30M");

        regionTask.processRecord(read1);

        SAMRecord read2 = createSamRecord(
                READ_ID, CHR_1, 1010, readBases, "30M", CHR_1, 1000, true,
                false, null, false, "30M");

        regionTask.processRecord(read2);

        assertEquals(1, posEvidence1.RefSupport);
        assertEquals(1, posEvidence1.ReadDepth);

        assertEquals(1, posEvidence1.RefSupport);
        assertEquals(1, posEvidence1.ReadDepth);

        read1 = createSamRecord(
                READ_ID, CHR_1, 1010, readBases, "30M", CHR_1, 1021, false,
                false, null, true, "30M");

        regionTask.processRecord(read1);

        read2 = createSamRecord(
                READ_ID, CHR_1, 1021, readBases, "30M", CHR_1, 1010, true,
                false, null, false, "30M");

        regionTask.processRecord(read2);

        assertEquals(2, posEvidence1.RefSupport);
        assertEquals(2, posEvidence1.ReadDepth);

        assertEquals(2, posEvidence1.RefSupport);
        assertEquals(2, posEvidence1.ReadDepth);
    }
}
