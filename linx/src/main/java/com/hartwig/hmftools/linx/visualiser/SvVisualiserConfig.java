package com.hartwig.hmftools.linx.visualiser;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;

import java.io.File;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

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
    public final RefGenomeCoordinates RefGenomeCoords;

    public final int Threads;
    public final boolean Debug;

    public final boolean PlotReportableEvents;
    public final boolean RestrictClusterByGene;

    private static final String PLOT_OUT = "plot_out";
    private static final String DATA_OUT = "data_out";
    private static final String CIRCOS = "circos";
    private static final String DEBUG = "debug";
    private static final String VIS_FILE_DIRECTORY = "vis_file_dir";
    private static final String PLOT_REPORTABLE = "plot_reportable";
    public static final String RESTRICT_CLUSTERS_BY_GENE = "restrict_cluster_by_gene";
    public static final String PLOT_CLUSTER_GENES = "plot_cluster_genes";

    private static final String THREADS = "threads";
    private static final String INCLUDE_LINE_ELEMENTS = "include_line_elements";

    public static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(PLOT_OUT, true, "Plot output directory");
        options.addOption(DATA_OUT, true, "Data output directory");
        options.addOption(VIS_FILE_DIRECTORY, true,
                "Path to all Linx vis files, used instead of specifying them individually");

        options.addOption(CIRCOS, true, "Path to circos binary");
        EnsemblDataCache.addEnsemblDir(options);

        options.addOption(DEBUG, false, "Enabled debug mode");
        options.addOption(THREADS, true, "Number of threads to use");
        options.addOption(INCLUDE_LINE_ELEMENTS, false, "Include line elements in chromosome plots");
        options.addOption(PLOT_REPORTABLE, false, "Plot all clusters with a fusion, disruption, AMP or DEL");

        options.addOption(PLOT_CLUSTER_GENES, false,
                "Show all genes linked to SVs in a cluster (only applicable when clusterId is specified)");

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

        RefGenomeVersion refGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));
        RefGenomeCoords = refGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        Debug = cmd.hasOption(DEBUG);

        IncludeLineElements = cmd.hasOption(INCLUDE_LINE_ELEMENTS);
        PlotReportableEvents = cmd.hasOption(PLOT_REPORTABLE);
        RestrictClusterByGene = cmd.hasOption(RESTRICT_CLUSTERS_BY_GENE);
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
