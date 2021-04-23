package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultIntValue;
import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

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
public interface SomaticFitConfig
{
    String SOMATIC_VARIANTS = "somatic_vcf";
    String SOMATIC_MIN_PEAK = "somatic_min_peak";
    String SOMATIC_MIN_TOTAL = "somatic_min_variants";
    String SOMATIC_MIN_PURITY = "somatic_min_purity";
    String SOMATIC_MIN_PURITY_SPREAD = "somatic_min_purity_spread";
    String SOMATIC_PENALTY_WEIGHT = "somatic_penalty_weight";
    String HIGHLY_DIPLOID_PERCENTAGE = "highly_diploid_percentage";
    String FORCE_SOMATIC_FIT = "force_somatic_fit";

    public double SOMATIC_MIN_PURITY_DEFAULT = 0.17;
    public double SOMATIC_MIN_PURITY_SPREAD_DEFAULT = 0.15;
    public int SOMATIC_MIN_PEAK_DEFAULT = 10;
    public int SOMATIC_MIN_VARIANTS_DEFAULT = 10;
    double SOMATIC_PENALTY_WEIGHT_DEFAULT = 1;
    double HIGHLY_DIPLOID_PERCENTAGE_DEFAULT = 0.97;

    double MIN_SOMATIC_TOTAL_READ_COUNT_PROPORTION = 0.6;
    double MAX_SOMATIC_TOTAL_READ_COUNT_PROPORTION = 1.4;

    static void addOptions(@NotNull Options options)
    {
        options.addOption(SOMATIC_MIN_PEAK, true,
                "Minimum number of somatic variants to consider a peak. Default " + SOMATIC_MIN_PEAK_DEFAULT);
        options.addOption(SOMATIC_MIN_TOTAL, true,
                "Minimum number of somatic variants required to assist highly diploid fits. Default " + SOMATIC_MIN_VARIANTS_DEFAULT);
        options.addOption(SOMATIC_MIN_PURITY, true,
                "Somatic fit will not be used if both somatic and fitted purities are less than this value. Default "
                        + SOMATIC_MIN_PURITY_DEFAULT);
        options.addOption(SOMATIC_MIN_PURITY_SPREAD, true,
                "Minimum spread within candidate purities before somatics can be used. Default " + SOMATIC_MIN_PURITY_SPREAD_DEFAULT);
        options.addOption(SOMATIC_VARIANTS, true, "Optional location of somatic variant vcf to assist fitting in highly-diploid samples.");
        options.addOption(SOMATIC_PENALTY_WEIGHT, true,
                "Proportion of somatic deviation to include in fitted purity score. Default " + SOMATIC_PENALTY_WEIGHT_DEFAULT);
        options.addOption(HIGHLY_DIPLOID_PERCENTAGE, true,
                "Proportion of genome that must be diploid before using somatic fit. Default " + HIGHLY_DIPLOID_PERCENTAGE_DEFAULT);
        options.addOption(FORCE_SOMATIC_FIT, false, "Fit from somatic VAFs only");
    }

    Optional<File> file();

    default boolean enabled()
    {
        return file().isPresent();
    }

    default int minTotalSomaticVariantAlleleReadCount()
    {
        return 5000;
    }

    default int minTotalSvFragmentCount()
    {
        return 1000;
    }

    int minTotalVariants();

    int minPeakVariants();

    double minSomaticPurity();

    double minSomaticPuritySpread();

    double somaticPenaltyWeight();

    double highlyDiploidPercentage();

    default double clonalityMaxPloidy()
    {
        return 10;
    }

    default double clonalityBinWidth()
    {
        return 0.05;
    }

    int minSomaticTotalReadCount();

    int maxSomaticTotalReadCount();

    boolean forceSomaticFit();

    @NotNull
    static SomaticFitConfig createSomaticConfig(
            @NotNull CommandLine cmd, final CommonConfig commonConfig, @NotNull final AmberData amberData) throws ParseException
    {
        final Optional<File> file;
        if(cmd.hasOption(SOMATIC_VARIANTS) || !commonConfig.sampleDirectory().isEmpty())
        {
            final String somaticFilename = cmd.getOptionValue(SOMATIC_VARIANTS,
                    commonConfig.sampleDirectory() + commonConfig.tumorSample() + ".sage.somatic.filtered.vcf.gz");

            final File somaticFile = new File(somaticFilename);
            if(!somaticFile.exists())
            {
                throw new ParseException("Unable to read somatic variants from: " + somaticFilename);
            }
            file = Optional.of(somaticFile);
        }
        else
        {
            file = Optional.empty();
            PPL_LOGGER.info("No somatic vcf supplied");
        }

        return ImmutableSomaticFitConfig.builder()
                .file(file)
                .minSomaticTotalReadCount((int) Math.floor(MIN_SOMATIC_TOTAL_READ_COUNT_PROPORTION * amberData.averageTumorDepth()))
                .maxSomaticTotalReadCount((int) Math.ceil(MAX_SOMATIC_TOTAL_READ_COUNT_PROPORTION * amberData.averageTumorDepth()))
                .minTotalVariants(defaultIntValue(cmd, SOMATIC_MIN_TOTAL, SOMATIC_MIN_VARIANTS_DEFAULT))
                .minPeakVariants(defaultIntValue(cmd, SOMATIC_MIN_PEAK, SOMATIC_MIN_PEAK_DEFAULT))
                .minSomaticPurity(defaultValue(cmd, SOMATIC_MIN_PURITY, SOMATIC_MIN_PURITY_DEFAULT))
                .minSomaticPuritySpread(defaultValue(cmd, SOMATIC_MIN_PURITY_SPREAD, SOMATIC_MIN_PURITY_SPREAD_DEFAULT))
                .somaticPenaltyWeight(defaultValue(cmd, SOMATIC_PENALTY_WEIGHT, SOMATIC_PENALTY_WEIGHT_DEFAULT))
                .highlyDiploidPercentage(defaultValue(cmd, HIGHLY_DIPLOID_PERCENTAGE, HIGHLY_DIPLOID_PERCENTAGE_DEFAULT))
                .forceSomaticFit(cmd.hasOption(FORCE_SOMATIC_FIT))
                .build();
    }
}
