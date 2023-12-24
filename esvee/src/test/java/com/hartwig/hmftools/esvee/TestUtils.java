package com.hartwig.hmftools.esvee;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.esvee.config.HMFConfig;
import com.hartwig.hmftools.esvee.models.Record;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class TestUtils
{
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

    public static Record createSAMRecord(final String sequence)
    {
        return createSAMRecord(sequence, 1);
    }

    public static Record createSAMRecord(final String sequence, final int position)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReferenceName("1");
        record.setReadName(String.valueOf(System.nanoTime()));
        record.setCigarString(sequence.length() + "M");
        record.setAlignmentStart(position);
        record.setReadString(sequence);
        record.setBaseQualityString("F".repeat(sequence.length()));
        record.setReadPairedFlag(true);
        record.setMateUnmappedFlag(true);
        return new Record(record);
    }
}
