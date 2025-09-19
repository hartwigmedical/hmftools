package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_MISMATCH_QUAL;
import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.SBX_DUPLEX_QUAL;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.READ_ID_GEN;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.redux.jitter.ConsensusMarker;
import com.hartwig.hmftools.redux.jitter.MicrosatelliteRead;
import com.hartwig.hmftools.redux.jitter.MicrosatelliteSite;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MsiJitterReadTest
{
    @Test
    public void testMicrosatelliteRead()
    {
        String repeat = "AAAA";
        MicrosatelliteSite msiSite = new MicrosatelliteSite(new ChrBaseRegion(CHR_1, 50, 53), repeat.getBytes());

        MicrosatelliteRead microsatelliteRead = new MicrosatelliteRead();

        String flankRefBases = "ACGTACCCCCGGGGGACGTA";
        String readBases = flankRefBases + repeat + flankRefBases;
        String cigar = format("%dM", readBases.length());

        SAMRecord read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 31, readBases, cigar, CHR_1, 200, false, false,
                null);

        ConsensusMarker consensusMarker = new ConsensusMarker.StandardConsensusMarker();

        microsatelliteRead.analyse(msiSite, read, consensusMarker);

        assertTrue(microsatelliteRead.isValidRead());

        // a read too close to the flanks on either side is invalid
        read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 46, readBases, cigar, CHR_1, 200, false, false,
                null);

        microsatelliteRead.analyse(msiSite, read, consensusMarker);

        assertFalse(microsatelliteRead.isValidRead());

        read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 8, readBases, cigar, CHR_1, 200, false, false,
                null);

        microsatelliteRead.analyse(msiSite, read, consensusMarker);

        assertFalse(microsatelliteRead.isValidRead());

        // low or medium quals within the repeat bounds is invalid
        read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 31, readBases, cigar, CHR_1, 200, false, false,
                null);

        read.getBaseQualities()[20] = SBX_DUPLEX_MISMATCH_QUAL;

        microsatelliteRead.analyse(msiSite, read, consensusMarker);

        assertFalse(microsatelliteRead.isValidRead());

        read.getBaseQualities()[20] = SBX_DUPLEX_QUAL;
        read.getBaseQualities()[23] = SBX_DUPLEX_MISMATCH_QUAL;

        microsatelliteRead.analyse(msiSite, read, consensusMarker);

        assertFalse(microsatelliteRead.isValidRead());

        // deletes covering the repeat
        cigar = "18M2D26M";
        read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 31, readBases, cigar, CHR_1, 200, false, false,
                null);

        microsatelliteRead.analyse(msiSite, read, consensusMarker);

        assertFalse(microsatelliteRead.isValidRead());

        // invalidating inserts
        cigar = "18M2I24M";
        read = createSamRecord(
                READ_ID_GEN.nextId(), CHR_1, 31, readBases, cigar, CHR_1, 200, false, false,
                null);

        microsatelliteRead.analyse(msiSite, read, consensusMarker);

        assertFalse(microsatelliteRead.isValidRead());

    }
}
