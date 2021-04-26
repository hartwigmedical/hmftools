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
public interface GermlineConfig
{

    Logger LOGGER = LogManager.getLogger(GermlineConfig.class);

    public static String GERMLINE_VARIANTS = "germline_vcf";

    static void addOptions(@NotNull Options options)
    {
        options.addOption(GERMLINE_VARIANTS, true, "Optional location of germline variants to enrich and process in driver catalog.");
    }

    Optional<File> file();

    default boolean enabled()
    {
        return file().isPresent();
    }

    @NotNull
    static GermlineConfig createGermlineConfig(@NotNull CommandLine cmd, final CommonConfig commonConfig) throws ParseException
    {
        final Optional<File> file;
        if(cmd.hasOption(GERMLINE_VARIANTS) || !commonConfig.sampleDirectory().isEmpty())
        {
            final String somaticFilename = cmd.getOptionValue(GERMLINE_VARIANTS,
                    commonConfig.sampleDirectory() + commonConfig.tumorSample() + ".sage.germline.filtered.vcf.gz");

            final File somaticFile = new File(somaticFilename);
            if(!somaticFile.exists())
            {
                throw new ParseException("Unable to read germline variants from: " + somaticFilename);
            }
            file = Optional.of(somaticFile);
        }
        else
        {
            file = Optional.empty();
            LOGGER.info("No germline vcf supplied");
        }

        return ImmutableGermlineConfig.builder().file(file).build();
    }
}
