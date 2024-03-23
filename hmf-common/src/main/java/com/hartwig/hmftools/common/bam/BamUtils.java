package com.hartwig.hmftools.common.bam;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static htsjdk.samtools.ValidationStringency.DEFAULT_STRINGENCY;

import java.io.File;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public final class BamUtils
{
    public static final String BAM_VALIDATION_STRINGENCY = "bam_validation";
    public static final String BAM_VALIDATION_STRINGENCY_DESC = "BAM validation: STRICT (default), LENIENT or SILENT";

    public static void addValidationStringencyOption(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(BAM_VALIDATION_STRINGENCY, false, BAM_VALIDATION_STRINGENCY_DESC, DEFAULT_STRINGENCY.toString());
    }

    public static ValidationStringency validationStringency(final ConfigBuilder configBuilder)
    {
        return ValidationStringency.valueOf(configBuilder.getValue(BAM_VALIDATION_STRINGENCY));
    }

    public static RefGenomeVersion deriveRefGenomeVersion(final String bamFile)
    {
        // assumes file exist and has a valid header
        SamReader samReader = SamReaderFactory.makeDefault().open(new File(bamFile));
        String firstChromosome = samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceName();
        return firstChromosome.startsWith(CHR_PREFIX) ? V38 : V37;
    }

    public static RefGenomeVersion deriveRefGenomeVersion(final SamReader samReader)
    {
        // assumes file exist and has a valid header
        String firstChromosome = samReader.getFileHeader().getSequenceDictionary().getSequences().get(0).getSequenceName();
        return firstChromosome.startsWith(CHR_PREFIX) ? V38 : V37;
    }
}
