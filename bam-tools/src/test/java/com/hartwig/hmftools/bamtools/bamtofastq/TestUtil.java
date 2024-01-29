package com.hartwig.hmftools.bamtools.bamtofastq;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

import java.util.List;
import java.util.function.BiConsumer;

import com.google.common.collect.Lists;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMRecord;

public abstract class TestUtil
{
    private static final int READ_LENGTH = 100;
    private static final String MATE_CHR = CHR_2;
    private static final String READ_BASES = "A".repeat(READ_LENGTH);
    private static final String CIGAR = READ_LENGTH + "M";
    private static final String CONSENSUS_ATTRIBUTE_VALUE = "10;10";

    public static SAMRecord unmappedSamRecord(final String readName)
    {
        return createSamRecord(readName, NO_CHROMOSOME_NAME, NO_POSITION, READ_BASES, NO_CIGAR, NO_CHROMOSOME_NAME, NO_POSITION, false, false, null);
    }

    public static SAMRecord mappedPrimarySamRecord(final String readName, final String chromosome, int readStart)
    {
        return createSamRecord(readName, chromosome, readStart, READ_BASES, CIGAR, MATE_CHR, 1, false, false, null);
    }

    public static SAMRecord supplementarySamRecord(final String readName, final String chromosome, int readStart)
    {
        return createSamRecord(readName, chromosome, readStart, READ_BASES, CIGAR, MATE_CHR, 1, false, true, null);
    }

    public static SAMRecord consensusSamRecord(final String readName, final String chromosome, int readStart)
    {
        SAMRecord samRecord = mappedPrimarySamRecord(readName, chromosome, readStart);
        samRecord.setAttribute(CONSENSUS_READ_ATTRIBUTE, CONSENSUS_ATTRIBUTE_VALUE);
        return samRecord;
    }

    public static SAMRecord duplicateSamRecord(final String readName, final String chromosome, int readStart)
    {
        SAMRecord samRecord = mappedPrimarySamRecord(readName, chromosome, readStart);
        samRecord.setDuplicateReadFlag(true);
        return samRecord;
    }

    public static Pair<SAMRecord, SAMRecord> pairedMappedPrimarySamRecords(final String readName, final String chromosome1, int readStart1,
            final String chromosome2, int readStart2)
    {
        SAMRecord read1 =
                createSamRecord(readName, chromosome1, readStart1, READ_BASES, CIGAR, chromosome2, readStart2, false, false, null);
        read1.setFirstOfPairFlag(true);
        read1.setSecondOfPairFlag(false);

        SAMRecord read2 =
                createSamRecord(readName, chromosome2, readStart2, READ_BASES, CIGAR, chromosome1, readStart1, false, false, null);
        read2.setFirstOfPairFlag(false);
        read2.setSecondOfPairFlag(true);

        return Pair.of(read1, read2);
    }

    public static class ReadPairCollector implements BiConsumer<SAMRecord, SAMRecord>
    {
        public final List<Pair<SAMRecord, SAMRecord>> ReadPairs;

        public ReadPairCollector()
        {
            ReadPairs = Lists.newArrayList();
        }

        @Override
        public void accept(final SAMRecord read1, final SAMRecord read2)
        {
            ReadPairs.add(Pair.of(read1, read2));
        }
    }
}
