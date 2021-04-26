package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.util.Optional;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ChartConfig
{
    String CIRCOS = "circos";
    String DISABLE = "no_charts";

    static void addOptions(@NotNull Options options)
    {
        options.addOption(CIRCOS, true, "Location of circos binary");
        options.addOption(DISABLE, false, "Disable charts");
    }

    boolean enabled();

    Optional<String> circosBinary();

    String plotDirectory();

    String circosDirectory();

    @NotNull
    static ChartConfig createCircosConfig(@NotNull CommandLine cmd, @NotNull CommonConfig config)
    {
        return ImmutableChartConfig.builder()
                .enabled(!cmd.hasOption(DISABLE))
                .plotDirectory(config.outputDirectory() + File.separator + "plot")
                .circosDirectory(config.outputDirectory() + File.separator + "circos")
                .circosBinary(cmd.hasOption(CIRCOS) ? Optional.of(cmd.getOptionValue(CIRCOS)) : Optional.empty())
                .build();
    }
}
