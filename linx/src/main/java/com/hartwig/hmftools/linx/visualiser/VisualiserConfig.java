package com.hartwig.hmftools.linx.visualiser;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.SpecificRegions.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.region.SpecificRegions.SPECIFIC_REGIONS_DESC;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.linx.LinxConfig.GERMLINE;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.linx.visualiser.circos.Highlights;

import org.apache.commons.cli.ParseException;

public class VisualiserConfig
{
    public final String Sample;
    public final String SampleDataDir;

    public final String AmberDir;
    public final String CobaltDir;

    public final boolean UseCohortFiles;
    public final boolean IsGermline;

    public final String OutputConfPath;
    public final String OutputPlotPath;
    public final String CircosBin;
    public final String EnsemblDataDir;
    public final RefGenomeVersion RefGenVersion;

    public final RefGenomeCoordinates RefGenomeCoords;

    public final int Threads;
    public final boolean Debug;

    // filters and plotting options
    public final boolean IncludeFragileSites;
    public final boolean IncludeLineElements;
    public final List<Integer> ClusterIds;
    public final List<Integer> ChainIds;
    public final List<String> Chromosomes;
    public final Set<String> Genes;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean PlotReportableEvents;
    public final boolean PlotClusterGenes;
    public final boolean RestrictClusterByGene;

    private static final String VIS_FILE_DIRECTORY = "vis_file_dir";
    private static final String AMBER_DIRECTORY = "amber_dir";
    private static final String COBALT_DIRECTORY = "cobalt_dir";
    private static final String LOAD_COHORT_FILES = "load_cohort_files";
    private static final String CLUSTER_IDS = "clusterId";
    private static final String CHAIN_IDS = "chainId";
    private static final String CHROMOSOMES = "chromosome";
    private static final String GENE = "gene";
    private static final String PLOT_OUT = "plot_out";
    private static final String DATA_OUT = "data_out";
    private static final String CIRCOS = "circos";
    private static final String DEBUG = "debug";
    private static final String PLOT_REPORTABLE = "plot_reportable";
    public static final String RESTRICT_CLUSTERS_BY_GENE = "restrict_cluster_by_gene";
    public static final String PLOT_CLUSTER_GENES = "plot_cluster_genes";

    private static final String INCLUDE_FRAGILE_SITES = "include_fragile_sites";
    private static final String INCLUDE_LINE_ELEMENTS = "include_line_elements";

    private static final String DELIM = ",";

    public VisualiserConfig(final ConfigBuilder configBuilder) throws IOException, ParseException
    {
        Sample = configBuilder.getValue(SAMPLE);
        SampleDataDir = checkAddDirSeparator(configBuilder.getValue(VIS_FILE_DIRECTORY));
        CobaltDir = checkAddDirSeparator(configBuilder.getValue(COBALT_DIRECTORY));
        AmberDir = checkAddDirSeparator(configBuilder.getValue(AMBER_DIRECTORY));

        OutputPlotPath = checkAddDirSeparator(configBuilder.getValue(PLOT_OUT, SampleDataDir + "plot/"));
        OutputConfPath = checkAddDirSeparator(configBuilder.getValue(DATA_OUT, SampleDataDir + "data/"));
        CircosBin = configBuilder.getValue(CIRCOS);
        UseCohortFiles = configBuilder.hasFlag(LOAD_COHORT_FILES);
        IsGermline = configBuilder.hasFlag(GERMLINE);

        RefGenVersion = RefGenomeVersion.from(configBuilder);
        EnsemblDataDir = configBuilder.getValue(ENSEMBL_DATA_DIR);

        SpecificRegions = loadSpecificRegions(configBuilder);
        ClusterIds = parseClusterIds(configBuilder);
        ChainIds = parseChainIds(configBuilder);
        Chromosomes = parseChromosomes(configBuilder);
        Genes = Sets.newHashSet();

        if(configBuilder.hasValue(GENE))
        {
            String geneStr = configBuilder.getValue(GENE);
            String delim = geneStr.contains(ITEM_DELIM) ? ITEM_DELIM : DELIM;
            Arrays.stream(geneStr.split(delim)).forEach(x -> Genes.add(x));
        }

        File outputDir = new File(OutputPlotPath);
        if(!outputDir.exists() && !outputDir.mkdirs())
        {
            throw new IOException("Unable to write to plot directory " + OutputPlotPath);
        }

        outputDir = new File(OutputConfPath);
        if(!outputDir.exists() && !outputDir.mkdirs())
        {
            throw new IOException("Unable to write to data directory " + OutputConfPath);
        }

        RefGenomeCoords = RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        Highlights.populateKnownSites(RefGenVersion);

        Threads = parseThreads(configBuilder);
        Debug = configBuilder.hasFlag(DEBUG);

        IncludeFragileSites = configBuilder.hasFlag(INCLUDE_FRAGILE_SITES);
        IncludeLineElements = configBuilder.hasFlag(INCLUDE_LINE_ELEMENTS);
        PlotReportableEvents = configBuilder.hasFlag(PLOT_REPORTABLE);
        RestrictClusterByGene = configBuilder.hasFlag(RESTRICT_CLUSTERS_BY_GENE);
        PlotClusterGenes = configBuilder.hasFlag(PLOT_CLUSTER_GENES);
    }

    private static List<Integer> parseClusterIds(final ConfigBuilder configBuilder)
    {
        List<Integer> result = Lists.newArrayList();

        if(configBuilder.hasValue(CLUSTER_IDS))
        {
            Arrays.stream(configBuilder.getValue(CLUSTER_IDS).split(DELIM, -1)).forEach(x -> result.add(Integer.parseInt(x)));
        }

        return result;
    }

    private static List<Integer> parseChainIds(final ConfigBuilder configBuilder)
    {
        List<Integer> result = Lists.newArrayList();

        if(configBuilder.hasValue(CHAIN_IDS))
        {
            Arrays.stream(configBuilder.getValue(CHAIN_IDS).split(DELIM, -1)).forEach(x -> result.add(Integer.parseInt(x)));
        }

        return result;
    }

    private static List<String> parseChromosomes(final ConfigBuilder configBuilder)
    {
        List<String> result = Lists.newArrayList();
        if(configBuilder.hasValue(CHROMOSOMES))
        {
            final String contigs = configBuilder.getValue(CHROMOSOMES);

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

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);

        configBuilder.addPath(
                VIS_FILE_DIRECTORY, true, "Path to all Linx vis files, used instead of specifying them individually");

        configBuilder.addPath(AMBER_DIRECTORY, false, "Path to directory containing AMBER output");
        configBuilder.addPath(COBALT_DIRECTORY, false, "Path to directory containing COBALT output");

        configBuilder.addFlag(LOAD_COHORT_FILES, "Load Linx cohort rather than per-sample vis files");
        configBuilder.addFlag(GERMLINE, "Load Linx germline VIS files");
        addRefGenomeVersion(configBuilder);

        configBuilder.addConfigItem(PLOT_OUT, "Plot output directory, default is 'plot' in sample files directory");
        configBuilder.addConfigItem(DATA_OUT, "Data output directory, default is 'data' in sample files directory");

        configBuilder.addPath(CIRCOS, true, "Path to Circos binary");
        EnsemblDataCache.addEnsemblDir(configBuilder);

        // filters
        configBuilder.addConfigItem(GENE, "Show canonical transcript for genes (separated by ','");
        configBuilder.addConfigItem(CLUSTER_IDS, "Limit to specified cluster IDs (separated by ',')");
        configBuilder.addConfigItem(CHAIN_IDS, "Limit to specified chain IDs (separated by ','), requires clusterId to be specified");
        configBuilder.addConfigItem(CHROMOSOMES, "Limit to specified chromosomes (separated by ',')");
        configBuilder.addConfigItem(SPECIFIC_REGIONS, SPECIFIC_REGIONS_DESC);

        // options
        configBuilder.addFlag(RESTRICT_CLUSTERS_BY_GENE, "Only plot clusters with breakends in configured 'gene' list");
        configBuilder.addFlag(INCLUDE_FRAGILE_SITES, "Include fragile sites in chromosome plots");
        configBuilder.addFlag(INCLUDE_LINE_ELEMENTS, "Include line elements in chromosome plots");
        configBuilder.addFlag(PLOT_REPORTABLE, "Plot all clusters with a fusion, disruption, AMP or DEL");

        configBuilder.addFlag(PLOT_CLUSTER_GENES,
                "Show all genes linked to SVs in a cluster (only applicable when clusterId is specified)");

        // debug and threads
        configBuilder.addFlag(DEBUG, "Enabled debug mode");
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
