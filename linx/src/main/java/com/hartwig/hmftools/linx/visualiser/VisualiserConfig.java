package com.hartwig.hmftools.linx.visualiser;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS_DESC;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.common.utils.FileReaderUtils;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.core.util.FileUtils;
import org.jetbrains.annotations.NotNull;

public class VisualiserConfig
{
    public final String Sample;
    public final String SampleDataDir;
    public final boolean UseCohortFiles;

    public final String OutputConfPath;
    public final String OutputPlotPath;
    public final String CircosBin;
    public final String EnsemblDataDir;
    public final RefGenomeVersion RefGenVersion;

    public final RefGenomeCoordinates RefGenomeCoords;

    public final int Threads;
    public final boolean Debug;

    // filters and plotting options
    public final boolean IncludeLineElements;
    public final List<Integer> Clusters;
    public final List<String> Chromosomes;
    public final Set<String> Genes;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean PlotReportableEvents;
    public final boolean PlotClusterGenes;
    public final boolean RestrictClusterByGene;

    private static final String SAMPLE = "sample";
    private static final String VIS_FILE_DIRECTORY = "vis_file_dir";
    private static final String LOAD_COHORT_FILES = "load_cohort_files";
    private static final String CLUSTERS = "clusterId";
    private static final String CHROMOSOMES = "chromosome";
    private static final String GENE = "gene";
    private static final String PLOT_OUT = "plot_out";
    private static final String DATA_OUT = "data_out";
    private static final String CIRCOS = "circos";
    private static final String DEBUG = "debug";
    private static final String PLOT_REPORTABLE = "plot_reportable";
    public static final String RESTRICT_CLUSTERS_BY_GENE = "restrict_cluster_by_gene";
    public static final String PLOT_CLUSTER_GENES = "plot_cluster_genes";

    private static final String INCLUDE_LINE_ELEMENTS = "include_line_elements";

    private static final String DELIM = ",";
    private static final String ITEM_DELIM = ";";

    public VisualiserConfig(final CommandLine cmd) throws ParseException, IOException
    {
        final StringJoiner missingJoiner = new StringJoiner(", ");

        Sample = parameter(cmd, SAMPLE, missingJoiner);
        SampleDataDir = checkAddDirSeparator(cmd.getOptionValue(VIS_FILE_DIRECTORY));

        OutputPlotPath = cmd.hasOption(PLOT_OUT) ? cmd.getOptionValue(PLOT_OUT) : SampleDataDir + "plot/";
        OutputConfPath = cmd.hasOption(DATA_OUT) ? cmd.getOptionValue(DATA_OUT) : SampleDataDir + "data/";
        CircosBin = parameter(cmd, CIRCOS, missingJoiner);
        UseCohortFiles = cmd.hasOption(LOAD_COHORT_FILES);

        RefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));
        EnsemblDataDir = cmd.getOptionValue(ENSEMBL_DATA_DIR);

        SpecificRegions = ChrBaseRegion.loadSpecificRegions(cmd);
        Clusters = parseClusters(cmd);
        Chromosomes = parseChromosomes(cmd);
        Genes = Sets.newHashSet();

        if(cmd.hasOption(GENE))
        {
            String geneStr = cmd.getOptionValue(GENE);
            String delim = geneStr.contains(ITEM_DELIM) ? ITEM_DELIM : DELIM;
            Arrays.stream(geneStr.split(delim)).forEach(x -> Genes.add(x));
        }

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

        Threads = parseThreads(cmd);
        Debug = cmd.hasOption(DEBUG);

        IncludeLineElements = cmd.hasOption(INCLUDE_LINE_ELEMENTS);
        PlotReportableEvents = cmd.hasOption(PLOT_REPORTABLE);
        RestrictClusterByGene = cmd.hasOption(RESTRICT_CLUSTERS_BY_GENE);
        PlotClusterGenes = cmd.hasOption(PLOT_CLUSTER_GENES);
    }

    private static List<Integer> parseClusters(final CommandLine cmd) throws ParseException
    {
        List<Integer> result = Lists.newArrayList();

        if(cmd.hasOption(CLUSTERS))
        {
            Arrays.stream(cmd.getOptionValue(CLUSTERS).split(DELIM, -1)).forEach(x -> result.add(Integer.parseInt(x)));
        }

        return result;
    }

    private static List<String> parseChromosomes(final CommandLine cmd)
    {
        List<String> result = Lists.newArrayList();
        if(cmd.hasOption(CHROMOSOMES))
        {
            final String contigs = cmd.getOptionValue(CHROMOSOMES);

            if(contigs.equals("All"))
            {
                Arrays.stream(HumanChromosome.values()).forEach(x -> result.add(x.toString()));
            }
            else
            {
                Collections.addAll(result, contigs.split(","));
            }
        }
        return result;
    }


    public static Options createOptions()
    {
        final Options options = new Options();

        options.addOption(SAMPLE, true, "Sample name");
        options.addOption(VIS_FILE_DIRECTORY, true, "Path to all Linx vis files, used instead of specifying them individually");
        options.addOption(LOAD_COHORT_FILES, false, "Load Linx cohort rather than per-sample vis files");
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);

        options.addOption(PLOT_OUT, true, "Plot output directory, default is 'plot' in sample files directory");
        options.addOption(DATA_OUT, true, "Data output directory, default is 'data' in sample files directory");

        options.addOption(CIRCOS, true, "Path to Circos binary");
        EnsemblDataCache.addEnsemblDir(options);

        options.addOption(DEBUG, false, "Enabled debug mode");
        addThreadOptions(options);

        // filters
        options.addOption(GENE, true, "Show canonical transcript for genes (separated by ','");
        options.addOption(CLUSTERS, true, "Only generate image for specified comma separated clusters");
        options.addOption(CHROMOSOMES, true, "Only generate image for specified comma separated chromosomes");
        options.addOption(SPECIFIC_REGIONS, true, SPECIFIC_REGIONS_DESC);

        // options
        options.addOption(RESTRICT_CLUSTERS_BY_GENE, false, "Only plot clusters with breakends in configured 'gene' list");
        addEnsemblDir(options);
        options.addOption(INCLUDE_LINE_ELEMENTS, false, "Include line elements in chromosome plots");
        options.addOption(PLOT_REPORTABLE, false, "Plot all clusters with a fusion, disruption, AMP or DEL");

        options.addOption(PLOT_CLUSTER_GENES, false,
                "Show all genes linked to SVs in a cluster (only applicable when clusterId is specified)");

        return options;
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
