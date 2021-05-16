package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.linx.LinxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableFusion;
import com.hartwig.hmftools.linx.visualiser.data.VisSvData;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.Segment;
import com.hartwig.hmftools.linx.visualiser.data.VisCopyNumbers;
import com.hartwig.hmftools.linx.visualiser.data.VisExons;
import com.hartwig.hmftools.linx.visualiser.data.VisLinks;
import com.hartwig.hmftools.linx.visualiser.data.VisProteinDomains;
import com.hartwig.hmftools.linx.visualiser.data.VisSegments;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumberFile;
import com.hartwig.hmftools.linx.visualiser.file.VisFusionFile;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExonFile;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSegmentFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class SvVisualiserConfig
{
    public final String OutputConfPath;
    public final String OutputPlotPath;
    public final String CircosBin;

    public final boolean IncludeLineElements;

    public final int Threads;
    public final boolean Debug;

    public final boolean PlotReportableEvents;

    private static final String PLOT_OUT = "plot_out";
    private static final String DATA_OUT = "data_out";
    private static final String CIRCOS = "circos";
    private static final String DEBUG = "debug";
    private static final String VIS_FILE_DIRECTORY = "vis_file_dir";
    private static final String PLOT_REPORTABLE = "plot_reportable";

    private static final String THREADS = "threads";
    private static final String INCLUDE_LINE_ELEMENTS = "include_line_elements";

    public static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(PLOT_OUT, true, "Plot output directory");
        options.addOption(DATA_OUT, true, "Data output directory");
        options.addOption(VIS_FILE_DIRECTORY, true, "Path to all Linx vis files, used instead of specifying them individually");
        options.addOption(CIRCOS, true, "Path to circos binary");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache files");

        options.addOption(DEBUG, false, "Enabled debug mode");
        options.addOption(THREADS, true, "Number of threads to use");
        options.addOption(INCLUDE_LINE_ELEMENTS, false, "Include line elements in chromosome plots");
        options.addOption(PLOT_REPORTABLE, false, "Plot all clusters with a fusion, disruption, AMP or DEL");

        return options;
    }

    public SvVisualiserConfig(@NotNull final CommandLine cmd) throws ParseException, IOException
    {
        final StringJoiner missingJoiner = new StringJoiner(", ");

        OutputPlotPath = parameter(cmd, PLOT_OUT, missingJoiner);
        OutputConfPath = parameter(cmd, DATA_OUT, missingJoiner);
        CircosBin = parameter(cmd, CIRCOS, missingJoiner);

        final String missing = missingJoiner.toString();
        if (!missing.isEmpty())
        {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        File outputDir = new File(OutputPlotPath);
        if (!outputDir.exists() && !outputDir.mkdirs())
        {
            throw new IOException("Unable to write to plot directory " + OutputPlotPath);
        }

        outputDir = new File(OutputConfPath);
        if (!outputDir.exists() && !outputDir.mkdirs())
        {
            throw new IOException("Unable to write to data directory " + OutputConfPath);
        }

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        Debug = cmd.hasOption(DEBUG);

        IncludeLineElements = cmd.hasOption(INCLUDE_LINE_ELEMENTS);
        PlotReportableEvents = cmd.hasOption(PLOT_REPORTABLE);
    }

    public static String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing)
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
