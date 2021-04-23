package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.util.Optional;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface StructuralVariantConfig
{

    Logger LOGGER = LogManager.getLogger(StructuralVariantConfig.class);

    String STRUCTURAL_VARIANTS = "structural_vcf";
    String STRUCTURAL_VARIANT_RECOVERY = "sv_recovery_vcf";

    static void addOptions(@NotNull Options options)
    {
        options.addOption(STRUCTURAL_VARIANTS, true, "Optional location of structural variant vcf for more accurate segmentation.");
        options.addOption(STRUCTURAL_VARIANT_RECOVERY, true, "Optional location of failing structural variants that may be recovered.");
    }

    Optional<File> file();

    Optional<File> recoveryFile();

    @NotNull
    static StructuralVariantConfig createStructuralVariantConfig(
            @NotNull CommandLine cmd, @NotNull Options opt, final CommonConfig commonConfig) throws ParseException
    {
        final Optional<File> file;

        if(cmd.hasOption(STRUCTURAL_VARIANTS) || !commonConfig.sampleDirectory().isEmpty())
        {
            final String somaticFilename = cmd.getOptionValue(STRUCTURAL_VARIANTS,
                    commonConfig.sampleDirectory() + commonConfig.tumorSample() + ".gripss.somatic.filtered.vcf.gz");

            final File somaticFile = new File(somaticFilename);
            if(!somaticFile.exists())
            {
                throw new ParseException("Unable to read structural variants from: " + somaticFilename);
            }

            file = Optional.of(somaticFile);
        }
        else
        {
            file = Optional.empty();
            LOGGER.info("No structural vcf supplied");
        }

        final Optional<File> optionalRecoveryFile;

        if(cmd.hasOption(STRUCTURAL_VARIANT_RECOVERY) || !commonConfig.sampleDirectory().isEmpty())
        {
            final String somaticFilename = cmd.getOptionValue(STRUCTURAL_VARIANT_RECOVERY,
                    commonConfig.sampleDirectory() + commonConfig.tumorSample() + ".gripss.somatic.vcf.gz");

            final File recoveryFile = new File(somaticFilename);
            if(!recoveryFile.exists())
            {
                throw new ParseException("Unable to read structural variant recovery file from: " + somaticFilename);
            }
            optionalRecoveryFile = Optional.of(recoveryFile);
        }
        else
        {
            optionalRecoveryFile = Optional.empty();
        }

        return ImmutableStructuralVariantConfig.builder().file(file).recoveryFile(optionalRecoveryFile).build();
    }

}
