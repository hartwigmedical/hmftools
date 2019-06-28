package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlterations;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Exons;
import com.hartwig.hmftools.linx.visualiser.data.Link;
import com.hartwig.hmftools.linx.visualiser.data.Links;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomains;
import com.hartwig.hmftools.linx.visualiser.data.Segment;
import com.hartwig.hmftools.linx.visualiser.data.Segments;

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
    String PROTEIN_DOMAIN = "protein_domain";
    String LINK = "link";
    String CIRCOS = "circos";
    String THREADS = "threads";
    String DEBUG = "debug";
    String CLUSTERS = "clusterId";
    String CHROMOSOMES = "chromosome";
    String CNA = "cna";
    String EXON = "exon";
    String SCALE_EXON = "scale_exons";

    @NotNull
    String sample();

    @NotNull
    List<Segment> segments();

    @NotNull
    List<Link> links();

    @NotNull
    List<CopyNumberAlteration> copyNumberAlterations();

    @NotNull
    List<ProteinDomain> proteinDomain();

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

    boolean scaleExons();

    @NotNull
    List<Integer> clusters();

    @NotNull
    List<String> chromosomes();

    @NotNull
    static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(PLOT_OUT, true, "Plot output directory");
        options.addOption(DATA_OUT, true, "Data output directory");
        options.addOption(SAMPLE, true, "Sample name");
        options.addOption(SEGMENT, true, "Path to track file");
        options.addOption(LINK, true, "Path to link file");
        options.addOption(PROTEIN_DOMAIN, true, "Path to protein domain file");
        options.addOption(CIRCOS, true, "Path to circos binary");
        options.addOption(THREADS, true, "Number of threads to use");
        options.addOption(DEBUG, false, "Enabled debug mode");

        options.addOption(CLUSTERS, true, "Only generate image for specified comma separated clusters");
        options.addOption(CHROMOSOMES, true, "Only generate image for specified comma separated chromosomes");
        options.addOption(CNA, true, "Path to copy number alterations");
        options.addOption(EXON, true, "Path to exons");
        options.addOption(SCALE_EXON, false, "Scale exon positions instead of interpolating them");

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
        final String proteinDomainPath = parameter(cmd, PROTEIN_DOMAIN, missingJoiner);

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
        final List<ProteinDomain> proteinDomains =
                ProteinDomains.readProteinDomains(proteinDomainPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());

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
                .proteinDomain(proteinDomains)
                .copyNumberAlterations(cna)
                .circosBin(circos)
                .threads(Integer.valueOf(cmd.getOptionValue(THREADS, "1")))
                .debug(cmd.hasOption(DEBUG))
                .clusters(clusters(cmd))
                .chromosomes(chromosomes(cmd))
                .scaleExons(cmd.hasOption(SCALE_EXON))
                .build();
    }

    @NotNull
    static List<Integer> clusters(@NotNull final CommandLine cmd) throws ParseException
    {
        List<Integer> result = Lists.newArrayList();
        if (cmd.hasOption(CLUSTERS))
        {
            final String clusters = cmd.getOptionValue(CLUSTERS);
            for (String clusterId : clusters.split(","))
            {
                try
                {
                    result.add(Integer.valueOf(clusterId));
                } catch (NumberFormatException e)
                {
                    throw new ParseException(CLUSTERS + " should be comma separated integer values");
                }
            }

        }

        return result;
    }

    @NotNull
    static List<String> chromosomes(@NotNull final CommandLine cmd)
    {
        List<String> result = Lists.newArrayList();
        if (cmd.hasOption(CHROMOSOMES))
        {
            final String contigs = cmd.getOptionValue(CHROMOSOMES);
            Collections.addAll(result, contigs.split(","));
        }
        return result;
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
