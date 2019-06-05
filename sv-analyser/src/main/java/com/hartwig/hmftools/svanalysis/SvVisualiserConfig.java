package com.hartwig.hmftools.svanalysis;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.svanalysis.visualisation.data.CopyNumberAlteration;
import com.hartwig.hmftools.svanalysis.visualisation.data.CopyNumberAlterations;
import com.hartwig.hmftools.svanalysis.visualisation.data.Exon;
import com.hartwig.hmftools.svanalysis.visualisation.data.Exons;
import com.hartwig.hmftools.svanalysis.visualisation.data.Link;
import com.hartwig.hmftools.svanalysis.visualisation.data.Links;
import com.hartwig.hmftools.svanalysis.visualisation.data.Segment;
import com.hartwig.hmftools.svanalysis.visualisation.data.Segments;

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
public interface SvVisualiserConfig
{
    Logger LOGGER = LogManager.getLogger(SvVisualiserConfig.class);

    String PLOT_OUT = "plot_out";
    String DATA_OUT = "data_out";
    String SAMPLE = "sample";
    String SEGMENT = "segment";
    String LINK = "link";
    String CIRCOS = "circos";
    String THREADS = "threads";
    String DEBUG = "debug";
    String SINGLE_CLUSTER = "clusterId";
    String SINGLE_CHROMOSOME = "chromosome";
    String CNA = "cna";
    String EXON = "exon";

    @NotNull
    String sample();

    @NotNull
    List<Segment> segments();

    @NotNull
    List<Link> links();

    @NotNull
    List<CopyNumberAlteration> copyNumberAlterations();

    @NotNull
    List<Exon> exons();

    @NotNull
    String outputConfPath();

    @NotNull
    String outputPlotPath();

    @NotNull
    String circosBin();

    int threads();

    boolean debug();

    @Nullable
    Integer singleCluster();

    @Nullable
    String singleChromosome();

    @NotNull
    static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(PLOT_OUT, true, "Plot output directory");
        options.addOption(DATA_OUT, true, "Data output directory");
        options.addOption(SAMPLE, true, "Sample name");
        options.addOption(SEGMENT, true, "Path to track file");
        options.addOption(LINK, true, "Path to link file");
        options.addOption(CIRCOS, true, "Path to circos binary");
        options.addOption(THREADS, true, "Number of threads to use");
        options.addOption(DEBUG, false, "Enabled debug mode");

        options.addOption(SINGLE_CLUSTER, true, "Only generate image for single cluster");
        options.addOption(SINGLE_CHROMOSOME, true, "Only generate image for singe chromosome");
        options.addOption(CNA, true, "Location of copy number alterations (optional alternative to db)");
        options.addOption(EXON, true, "Location of exons");

        return options;
    }

    @NotNull
    static SvVisualiserConfig createConfig(@NotNull final CommandLine cmd) throws ParseException, IOException
    {
        final StringJoiner missingJoiner = new StringJoiner(", ");
        final String linkPath = parameter(cmd, LINK, missingJoiner);
        final String trackPath = parameter(cmd, SEGMENT, missingJoiner);
        final String cnaPath = parameter(cmd, CNA, missingJoiner);
        final String sample = parameter(cmd, SAMPLE, missingJoiner);
        final String plotOutputDir = parameter(cmd, PLOT_OUT, missingJoiner);
        final String dataOutputDir = parameter(cmd, DATA_OUT, missingJoiner);
        final String circos = parameter(cmd, CIRCOS, missingJoiner);
        final String exonPath = parameter(cmd, EXON, missingJoiner);

        final String missing = missingJoiner.toString();
        if (!missing.isEmpty())
        {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        final List<Link> links = Links.readLinks(linkPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<Exon> exons = Exons.readExons(exonPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<Segment> segments = Segments.readTracks(trackPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<CopyNumberAlteration> cna =
                CopyNumberAlterations.read(cnaPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());

        if (segments.isEmpty() && links.isEmpty())
        {
            LOGGER.warn("No structural variants found for sample {}", sample);
        }

        if (cna.isEmpty())
        {
            LOGGER.warn("No copy number alterations found for sample {}", sample);
        }

        File outputDir = new File(plotOutputDir);
        if (!outputDir.exists() && !outputDir.mkdirs())
        {
            throw new IOException("Unable to write to plot directory " + plotOutputDir);
        }

        outputDir = new File(dataOutputDir);
        if (!outputDir.exists() && !outputDir.mkdirs())
        {
            throw new IOException("Unable to write to data directory " + plotOutputDir);
        }

        return ImmutableSvVisualiserConfig.builder()
                .outputConfPath(dataOutputDir)
                .outputPlotPath(plotOutputDir)
                .segments(segments)
                .links(links)
                .sample(sample)
                .exons(exons)
                .copyNumberAlterations(cna)
                .circosBin(circos)
                .threads(Integer.valueOf(cmd.getOptionValue(THREADS, "1")))
                .debug(cmd.hasOption(DEBUG))
                .singleCluster(cmd.hasOption(SINGLE_CLUSTER) ? Integer.valueOf(cmd.getOptionValue(SINGLE_CLUSTER)) : null)
                .singleChromosome(cmd.hasOption(SINGLE_CHROMOSOME) ? cmd.getOptionValue(SINGLE_CHROMOSOME) : null)
                .build();
    }

    @NotNull
    static String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing)
    {
        final String value = cmd.getOptionValue(parameter);
        if (value == null)
        {
            missing.add(parameter);
            return "";
        }
        return value;
    }

}
