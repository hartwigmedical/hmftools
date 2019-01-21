package com.hartwig.hmftools.svanalysis;

import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.svanalysis.visualisation.Link;
import com.hartwig.hmftools.svanalysis.visualisation.LinkFile;
import com.hartwig.hmftools.svanalysis.visualisation.Track;
import com.hartwig.hmftools.svanalysis.visualisation.TrackFile;

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
public interface SvVisualiserConfig {

    Logger LOGGER = LogManager.getLogger(SvVisualiserConfig.class);

    String OUT_PATH = "out";
    String SAMPLE = "sample";
    String TRACK = "track";
    String LINK = "link";
    String CIRCOS = "circos";

    @NotNull
    String sample();

    @NotNull
    List<Track> tracks();

    @NotNull
    List<Link> links();

    @NotNull
    String outputConfPath();

    @NotNull
    String outputPlotPath();

    @NotNull
    String circosBin();


    default int maxTracks() {
        return tracks().stream().mapToInt(Track::track).max().orElse(0) + 1;
    }

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_PATH, true, "Output directory");
        options.addOption(SAMPLE, true, "Sample name");
        options.addOption(TRACK, true, "Path to track file");
        options.addOption(LINK, true, "Path to link file");
        options.addOption(CIRCOS, true, "Path to circos binary");

        return options;
    }

    @NotNull
    static SvVisualiserConfig createConfig(@NotNull final CommandLine cmd) throws ParseException, IOException {
        final StringJoiner missingJoiner = new StringJoiner(", ");
        final String linkPath = parameter(cmd, LINK, missingJoiner);
        final String trackPath = parameter(cmd, TRACK, missingJoiner);
        final String sample = parameter(cmd, SAMPLE, missingJoiner);
        final String outputDir = parameter(cmd, OUT_PATH, missingJoiner);
        final String circos = parameter(cmd, CIRCOS, missingJoiner);
        final String missing = missingJoiner.toString();

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        final List<Track> tracks = TrackFile.readTracks(trackPath);
        final List<Link> links = LinkFile.readLinks(linkPath);

        return ImmutableSvVisualiserConfig.builder()
                .outputConfPath(outputDir)
                .outputPlotPath(outputDir)
                .tracks(tracks)
                .links(links)
                .sample(sample)
                .circosBin(circos)
                .build();
    }

    @NotNull
    static String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing) {
        final String value = cmd.getOptionValue(parameter);
        if (value == null) {
            missing.add(parameter);
            return "";
        }
        return value;
    }
}
