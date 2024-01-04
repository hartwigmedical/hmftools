package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.read.Read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class TestUtils
{
    /*
    public static SVAConfig config()
    {
        return config(Map.of());
    }

    public static SVAConfig config(final Map<String, Object> overrideConfig)
    {
        final Map<String, Object> baseConfig = Map.of(
                "bam_file", "dummy",
                "ref_genome", "dummy",
                "ref_genome_index", "dummy",
                "junction_file", "dummy",
                "output_file", "dummy");

        final Map<String, Object> combinedConfig = new HashMap<>();
        combinedConfig.putAll(baseConfig);
        combinedConfig.putAll(overrideConfig);

        return null; // HMFConfig.load(combinedConfig, SVAConfig.class, ImmutableSVAConfig.builder());
    }
    */

    public static Read createSAMRecord(final String sequence)
    {
        return createSAMRecord(sequence, 1);
    }

    public static Read createSAMRecord(final String sequence, final int position)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReferenceName(CHR_1);
        record.setReadName(String.valueOf(System.nanoTime()));
        record.setCigarString(sequence.length() + "M");
        record.setAlignmentStart(position);
        record.setReadString(sequence);
        record.setBaseQualityString("F".repeat(sequence.length()));
        record.setReadPairedFlag(true);
        record.setMateUnmappedFlag(true);
        return new Read(record);
    }

    public static Read createSamRecord(final String readId, int readStart, final String readBases, final String cigar)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, readStart + 1000,
                false, false, null);
        return new Read(record);
    }
}