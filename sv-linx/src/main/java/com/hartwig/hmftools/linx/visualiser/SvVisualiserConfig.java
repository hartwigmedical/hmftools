package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.linx.LinxConfig.GENE_TRANSCRIPTS_DIR;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlterations;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Exons;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.Fusions;
import com.hartwig.hmftools.linx.visualiser.data.Link;
import com.hartwig.hmftools.linx.visualiser.data.Links;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomains;
import com.hartwig.hmftools.linx.visualiser.data.Segment;
import com.hartwig.hmftools.linx.visualiser.data.Segments;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumberFile;
import com.hartwig.hmftools.linx.visualiser.file.VisFusionFile;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExonFile;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSegmentFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile;

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
    String FUSION = "fusion";
    String LINK = "link";
    String CIRCOS = "circos";
    String DEBUG = "debug";
    String CLUSTERS = "clusterId";
    String CHROMOSOMES = "chromosome";
    String CNA = "cna";
    String EXON = "exon";
    String VIS_FILE_DIRECTORY = "vis_file_dir";

    String THREADS = "threads";
    String INCLUDE_LINE_ELEMENTS = "include_line_elements";
    String GENE = "gene";

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
    List<Fusion> fusions();

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

    boolean includeLineElements();

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
        options.addOption(SEGMENT, true, "Path to segment file - eg 'COLO829T.linx.vis_segments.tsv");
        options.addOption(LINK, true, "Path to sv-data file eg 'COLO829T.linx.vis_sv_data.tsv'");
        options.addOption(PROTEIN_DOMAIN, true, "Path to protein domain file - eg 'COLO829T.linx.vis_protein_domain.tsv'");
        options.addOption(FUSION, true, "Path to fusion file - eg 'COLO829T.linx.fusions.tsv'");
        options.addOption(CNA, true, "Path to copy number alteration file - eg 'COLO829T.linx.vis_copy_number.tsv'");
        options.addOption(EXON, true, "Path to exon file - eg 'COLO829T.linx.vis_gene_exon.tsv'");
        options.addOption(VIS_FILE_DIRECTORY, true, "Path to all Linx vis files, used instead of specifying them individually");
        options.addOption(CIRCOS, true, "Path to circos binary");

        options.addOption(GENE, true, "Add canonical transcriptions of supplied comma separated genes to image");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache files");
        options.addOption(CLUSTERS, true, "Only generate image for specified comma separated clusters");
        options.addOption(CHROMOSOMES, true, "Only generate image for specified comma separated chromosomes");

        options.addOption(DEBUG, false, "Enabled debug mode");
        options.addOption(THREADS, true, "Number of threads to use");
        options.addOption(INCLUDE_LINE_ELEMENTS, false, "Include line elements in chromosome plots");

        return options;
    }

    @NotNull
    static SvVisualiserConfig createConfig(@NotNull final CommandLine cmd) throws ParseException, IOException
    {
        final StringJoiner missingJoiner = new StringJoiner(", ");

        final String sample = parameter(cmd, SAMPLE, missingJoiner);
        final String visFilesDirectory = cmd.getOptionValue(VIS_FILE_DIRECTORY);

        final String linkPath = visFilesDirectory != null ?
                VisSvDataFile.generateFilename(visFilesDirectory, sample) : parameter(cmd, LINK, missingJoiner);

        final String trackPath = visFilesDirectory != null ?
                VisSegmentFile.generateFilename(visFilesDirectory, sample) : parameter(cmd, SEGMENT, missingJoiner);

        final String cnaPath = visFilesDirectory != null ?
                VisCopyNumberFile.generateFilename(visFilesDirectory, sample) : parameter(cmd, CNA, missingJoiner);

        final String exonPath = visFilesDirectory != null ?
                VisGeneExonFile.generateFilename(visFilesDirectory, sample) : parameter(cmd, EXON, missingJoiner);

        final String proteinDomainPath = visFilesDirectory != null ?
                VisProteinDomainFile.generateFilename(visFilesDirectory, sample) : parameter(cmd, PROTEIN_DOMAIN, missingJoiner);

        final String fusionPath = visFilesDirectory != null ?
                VisFusionFile.generateFilename(visFilesDirectory, sample) : parameter(cmd, FUSION, missingJoiner);

        final String plotOutputDir = parameter(cmd, PLOT_OUT, missingJoiner);
        final String dataOutputDir = parameter(cmd, DATA_OUT, missingJoiner);
        final String circos = parameter(cmd, CIRCOS, missingJoiner);

        final String missing = missingJoiner.toString();
        if (!missing.isEmpty())
        {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        final List<Fusion> fusions = Fusions.fromFile(fusionPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<Link> links = Links.readLinks(linkPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<Exon> exons = Exons.readExons(exonPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<Segment> segments = Segments.readTracks(trackPath).stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<CopyNumberAlteration> cna = CopyNumberAlterations.read(cnaPath)
                .stream().filter(x -> x.sampleId().equals(sample)).collect(toList());
        final List<ProteinDomain> proteinDomains = ProteinDomains.readProteinDomains(proteinDomainPath, fusions)
                .stream()
                .filter(x -> x.sampleId().equals(sample))
                .collect(toList());

        final List<Integer> clusterIds = clusters(cmd);
        exons.addAll(additionalExons(cmd, exons, clusterIds));

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
                .fusions(fusions)
                .copyNumberAlterations(cna)
                .circosBin(circos)
                .threads(Integer.valueOf(cmd.getOptionValue(THREADS, "1")))
                .debug(cmd.hasOption(DEBUG))
                .clusters(clusterIds)
                .chromosomes(chromosomes(cmd))
                .includeLineElements(cmd.hasOption(INCLUDE_LINE_ELEMENTS))
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

    @NotNull
    static List<Exon> additionalExons(@NotNull final CommandLine cmd, @NotNull final List<Exon> currentExons,
            @NotNull final List<Integer> clusterIds)
    {
        final List<Exon> exonList = Lists.newArrayList();

        if (!cmd.hasOption(GENE))
            return exonList;

        final String sampleId = cmd.getOptionValue(SAMPLE);
        final List<Integer> allClusterIds = clusterIds.isEmpty() ? Lists.newArrayList(0) : clusterIds;

        final String[] geneList = cmd.getOptionValue(GENE).split(",");

        EnsemblDataCache geneTransCache = null;
        Map<String, HmfTranscriptRegion> geneMap = null;

        if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
        {
            geneTransCache = new EnsemblDataCache(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR), RefGenomeVersion.HG19);
            geneTransCache.setRequiredData(true, false, false, true);
            geneTransCache.load(false);
        }
        else
        {
            geneMap = HmfGenePanelSupplier.allGenesMap37();
        }

        for (final String geneName : geneList)
        {
            if (currentExons.stream().anyMatch(x -> x.gene().equals(geneName)))
                continue;

            LOGGER.info("loading exon data for additional gene({}}", geneName);

            if(geneTransCache != null)
            {
                EnsemblGeneData geneData = geneTransCache.getGeneDataByName(geneName);
                TranscriptData transcriptData = geneData != null ? geneTransCache.getTranscriptData(geneData.GeneId, "") : null;

                if (transcriptData == null)
                {
                    LOGGER.warn("data not found for specified gene({})", geneName);
                    continue;
                }

                for (Integer clusterId : allClusterIds)
                {
                    exonList.addAll(Exons.extractExonList(sampleId, clusterId, geneData, transcriptData));
                }
            }
            else
            {
                HmfTranscriptRegion hmfGene = geneMap.get(geneName);
                if (hmfGene == null)
                {
                    LOGGER.warn("data not found for specified gene({})", geneName);
                }
                else
                {
                    for (Integer clusterId : allClusterIds)
                    {
                        exonList.addAll(Exons.extractExonList(sampleId, clusterId, hmfGene));
                    }
                }
            }
        }

        return exonList;
    }

}
