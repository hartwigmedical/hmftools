package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiserConfig.parameter;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableFusion;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.Segment;
import com.hartwig.hmftools.linx.visualiser.data.VisCopyNumbers;
import com.hartwig.hmftools.linx.visualiser.data.VisExons;
import com.hartwig.hmftools.linx.visualiser.data.VisLinks;
import com.hartwig.hmftools.linx.visualiser.data.VisProteinDomains;
import com.hartwig.hmftools.linx.visualiser.data.VisSegments;
import com.hartwig.hmftools.linx.visualiser.data.VisSvData;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumberFile;
import com.hartwig.hmftools.linx.visualiser.file.VisFusionFile;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExonFile;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSegmentFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class SampleData
{
    public final String Sample;

    public final List<Segment> Segments;
    public final List<VisSvData> SvData;
    public final List<CopyNumberAlteration> CopyNumberAlterations;
    public final List<ProteinDomain> ProteinDomains;
    public final List<Fusion> Fusions;
    public final List<Exon> Exons;
    public final List<Integer> Clusters;

    public final List<String> Chromosomes;

    private final String mSampleDataDir;

    private static final String SAMPLE = "sample";
    private static final String VIS_FILE_DIRECTORY = "vis_file_dir";

    private static final String SEGMENT = "segment";
    private static final String PROTEIN_DOMAIN = "protein_domain";
    private static final String FUSION = "fusion";
    private static final String LINK = "link";

    private static final String CLUSTERS = "clusterId";
    private static final String CHROMOSOMES = "chromosome";
    private static final String CNA = "cna";
    private static final String EXON = "exon";
    private static final String GENE = "gene";

    private static final String PLOT_CLUSTER_GENES = "plot_cluster_genes";

    public SampleData(final CommandLine cmd) throws ParseException, IOException
    {
        final StringJoiner missingJoiner = new StringJoiner(", ");

        Sample = parameter(cmd, SAMPLE, missingJoiner);
        mSampleDataDir = cmd.getOptionValue(VIS_FILE_DIRECTORY);

        final String svDataFile = mSampleDataDir != null ?
                VisSvDataFile.generateFilename(mSampleDataDir, Sample) : parameter(cmd, LINK, missingJoiner);

        final String linksFile = mSampleDataDir != null ?
                VisSegmentFile.generateFilename(mSampleDataDir, Sample) : parameter(cmd, SEGMENT, missingJoiner);

        final String cnaFile = mSampleDataDir != null ?
                VisCopyNumberFile.generateFilename(mSampleDataDir, Sample) : parameter(cmd, CNA, missingJoiner);

        final String geneExonFile = mSampleDataDir != null ?
                VisGeneExonFile.generateFilename(mSampleDataDir, Sample) : parameter(cmd, EXON, missingJoiner);

        final String proteinDomainFile = mSampleDataDir != null ?
                VisProteinDomainFile.generateFilename(mSampleDataDir, Sample) : parameter(cmd, PROTEIN_DOMAIN, missingJoiner);

        final String fusionFile = mSampleDataDir != null ?
                VisFusionFile.generateFilename(mSampleDataDir, Sample) : parameter(cmd, FUSION, missingJoiner);

        final String missing = missingJoiner.toString();
        if (!missing.isEmpty())
        {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        Fusions = loadFusions(fusionFile).stream().filter(x -> x.sampleId().equals(Sample)).collect(toList());
        SvData = VisLinks.readSvData(svDataFile).stream().filter(x -> x.sampleId().equals(Sample)).collect(toList());

        Exons = VisExons.readExons(geneExonFile).stream().filter(x -> x.sampleId().equals(Sample)).collect(toList());
        Segments = VisSegments.readTracks(linksFile).stream().filter(x -> x.sampleId().equals(Sample)).collect(toList());

        CopyNumberAlterations = VisCopyNumbers.read(cnaFile)
                .stream().filter(x -> x.sampleId().equals(Sample)).collect(toList());

        ProteinDomains = VisProteinDomains.readProteinDomains(proteinDomainFile, Fusions)
                .stream()
                .filter(x -> x.sampleId().equals(Sample))
                .collect(toList());

        Clusters = parseClusters(cmd);
        Chromosomes = parseChromosomes(cmd);

        final Set<String> geneList = Sets.newHashSet();

        if(cmd.hasOption(GENE))
        {
            Arrays.stream(cmd.getOptionValue(GENE).split(",")).forEach(x -> geneList.add(x));
        }

        if(cmd.hasOption(PLOT_CLUSTER_GENES) && !Clusters.isEmpty() && mSampleDataDir != null)
        {
            List<LinxSvAnnotation> svAnnotations = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(mSampleDataDir, Sample));

            for(LinxSvAnnotation svAnnotation : svAnnotations)
            {
                if(Clusters.contains(svAnnotation.clusterId()))
                {
                    if(!svAnnotation.geneStart().isEmpty())
                        Arrays.stream(svAnnotation.geneStart().split(ITEM_DELIM)).forEach(x -> geneList.add(x));

                    if(!svAnnotation.geneEnd().isEmpty())
                        Arrays.stream(svAnnotation.geneEnd().split(ITEM_DELIM)).forEach(x -> geneList.add(x));
                }
            }
        }

        Exons.addAll(additionalExons(geneList, cmd, Exons, Clusters));

        if(Segments.isEmpty() && SvData.isEmpty())
        {
            VIS_LOGGER.warn("No structural variants found for sample {}", Sample);
        }

        if(CopyNumberAlterations.isEmpty())
        {
            VIS_LOGGER.warn("No copy number alterations found for sample {}", Sample);
        }
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(SAMPLE, true, "Sample name");
        options.addOption(VIS_FILE_DIRECTORY, true, "Path to all Linx vis files, used instead of specifying them individually");
        options.addOption(SEGMENT, true, "Path to segment file - eg 'COLO829T.linx.vis_segments.tsv");
        options.addOption(LINK, true, "Path to sv-data file eg 'COLO829T.linx.vis_sv_data.tsv'");
        options.addOption(PROTEIN_DOMAIN, true, "Path to protein domain file - eg 'COLO829T.linx.vis_protein_domain.tsv'");
        options.addOption(FUSION, true, "Path to fusion file - eg 'COLO829T.linx.fusions.tsv'");
        options.addOption(CNA, true, "Path to copy number alteration file - eg 'COLO829T.linx.vis_copy_number.tsv'");
        options.addOption(EXON, true, "Path to exon file - eg 'COLO829T.linx.vis_gene_exon.tsv'");
        options.addOption(GENE, true, "Add canonical transcriptions of supplied comma separated genes to image");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache files");
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - accepts 37 (default), or 38");
        options.addOption(CLUSTERS, true, "Only generate image for specified comma separated clusters");
        options.addOption(CHROMOSOMES, true, "Only generate image for specified comma separated chromosomes");

        options.addOption(PLOT_CLUSTER_GENES, false,
                "Show all genes linked to SVs in a cluster (only applicable when clusterId is specified)");
    }

    public Set<Integer> findReportableClusters()
    {
        Set<Integer> clusterIds = Sets.newHashSet();

        Fusions.stream().forEach(x -> clusterIds.add(x.clusterId()));

        if(mSampleDataDir != null)
        {
            try
            {
                // reportable disruptions
                final List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(mSampleDataDir, Sample));

                final List<Integer> svIds = breakends.stream()
                        .filter(x -> x.reportedDisruption()).map(x -> x.svId()).collect(toList());

                for(Integer svId : svIds)
                {
                    VisSvData svData = SvData.stream().filter(x -> x.svId() == svId).findFirst().orElse(null);
                    if(svData != null)
                         clusterIds.add(svData.clusterId());
                }

                final List<LinxDriver> drivers = LinxDriver.read(LinxDriver.generateFilename(mSampleDataDir, Sample));
                drivers.stream().filter(x -> x.clusterId() >= 0).forEach(x -> clusterIds.add(x.clusterId()));
            }
            catch(Exception e)
            {
                VIS_LOGGER.error("sample({}) could not read breakends or drivers: {}", Sample, e.toString());
            }
        }

        return clusterIds;
    }

    private List<Fusion> loadFusions(final String fileName) throws IOException
    {
        if(!Files.exists(Paths.get(fileName)))
            return Lists.newArrayList();

        final List<VisFusionFile> visFusions = VisFusionFile.read(fileName);

        return visFusions.stream().map(x -> ImmutableFusion.builder()
                .sampleId(x.SampleId)
                .reportable(x.Reportable)
                .clusterId(x.ClusterId)
                .geneUp(x.GeneNameUp)
                .transcriptUp(x.TranscriptUp)
                .chromosomeUp(x.ChrUp)
                .positionUp(x.PosUp)
                .strandUp(x.StrandUp)
                .regionTypeUp(x.RegionTypeUp)
                .fusedExonUp(x.FusedExonUp)
                .geneDown(x.GeneNameDown)
                .transcriptDown(x.TranscriptDown)
                .chromosomeDown(x.ChrDown)
                .positionDown(x.PosDown)
                .strandDown(x.StrandDown)
                .regionTypeDown(x.RegionTypeDown)
                .fusedExonDown(x.FusedExonDown)
                .build()).collect(Collectors.toList());
    }

    private static List<Integer> parseClusters(final CommandLine cmd) throws ParseException
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

    private static List<String> parseChromosomes(final CommandLine cmd)
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

    private static List<Exon> additionalExons(
            final Set<String> geneList, final CommandLine cmd, final List<Exon> currentExons, final List<Integer> clusterIds)
    {
        final List<Exon> exonList = Lists.newArrayList();

        if(geneList.isEmpty())
            return exonList;

        final String sampleId = cmd.getOptionValue(SAMPLE);
        final List<Integer> allClusterIds = clusterIds.isEmpty() ? Lists.newArrayList(0) : clusterIds;

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));

        EnsemblDataCache geneTransCache = null;
        Map<String, HmfTranscriptRegion> geneMap = null;

        if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
        {
            geneTransCache = new EnsemblDataCache(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR), refGenomeVersion);
            geneTransCache.setRequiredData(true, false, false, true);
            geneTransCache.load(false);
        }
        else
        {
            geneMap = refGenomeVersion.is37() ? HmfGenePanelSupplier.allGenesMap37() : HmfGenePanelSupplier.allGenesMap38();
        }

        for (final String geneName : geneList)
        {
            if (currentExons.stream().anyMatch(x -> x.gene().equals(geneName)))
                continue;

            VIS_LOGGER.info("loading exon data for additional gene({})", geneName);

            if(geneTransCache != null)
            {
                EnsemblGeneData geneData = geneTransCache.getGeneDataByName(geneName);
                TranscriptData transcriptData = geneData != null ? geneTransCache.getTranscriptData(geneData.GeneId, "") : null;

                if (transcriptData == null)
                {
                    VIS_LOGGER.warn("data not found for specified gene({})", geneName);
                    continue;
                }

                for (Integer clusterId : allClusterIds)
                {
                    exonList.addAll(VisExons.extractExonList(sampleId, clusterId, geneData, transcriptData));
                }
            }
            else
            {
                HmfTranscriptRegion hmfGene = geneMap.get(geneName);
                if (hmfGene == null)
                {
                    VIS_LOGGER.warn("data not found for specified gene({})", geneName);
                }
                else
                {
                    for (Integer clusterId : allClusterIds)
                    {
                        exonList.addAll(VisExons.extractExonList(sampleId, clusterId, hmfGene));
                    }
                }
            }
        }

        return exonList;
    }

}
