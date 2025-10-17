package com.hartwig.hmftools.common.test;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.PHRED_OFFSET;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_T0_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_TP_TAG;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_MAP_QUAL;

import htsjdk.samtools.SAMRecord;

public class SeqTechTestUtils
{
    public static SAMRecord buildUltimaRead(
            final String readId, final String chromosome, final int readStart, final String readBases, final byte[] qualities,
            final byte[] tpValues, final byte[] t0Values)
    {
        String cigar = format("%dM", qualities.length);

        // SAMRecord record = buildSamRecord(readStart, cigar, readBases, qualities);
        final SAMRecord record = new SAMRecord(null);
        record.setReadName(readId);
        record.setReferenceName(CHR_1);
        record.setAlignmentStart(readStart);
        record.setCigarString(cigar);
        record.setReadString(readBases);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualities(qualities);
        record.setMappingQuality(DEFAULT_MAP_QUAL);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(false);
        record.setReadPairedFlag(false);
        // record.setInferredInsertSize(600);
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, 0);

        record.setAttribute(ULTIMA_TP_TAG, tpValues);

        for(int i = 0; i < t0Values.length; ++i)
        {
            t0Values[i] += PHRED_OFFSET;
        }

        record.setAttribute(ULTIMA_T0_TAG, new String(t0Values));
        return record;
    }

}
