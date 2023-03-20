package com.hartwig.hmftools.common.samtools;

import static htsjdk.samtools.ValidationStringency.DEFAULT_STRINGENCY;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import htsjdk.samtools.ValidationStringency;

public final class BamUtils
{
    public static final String BAM_VALIDATION_STRINGENCY = "bam_validation";
    public static final String BAM_VALIDATION_STRINGENCY_DESC = "BAM validation: STRICT (default), LENIENT or SILENT";

    public static void addValidationStringencyOption(final Options options)
    {
        options.addOption(BAM_VALIDATION_STRINGENCY, true, BAM_VALIDATION_STRINGENCY_DESC);
    }

    public static ValidationStringency validationStringency(final CommandLine cmd)
    {
        return ValidationStringency.valueOf(cmd.getOptionValue(BAM_VALIDATION_STRINGENCY, DEFAULT_STRINGENCY.toString()));
    }

}
