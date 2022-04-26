package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS_DESC;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiserConfig.PLOT_CLUSTER_GENES;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiserConfig.RESTRICT_CLUSTERS_BY_GENE;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiserConfig.parameter;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_COPY_NUMBER_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_FUSIONS_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_GENE_EXONS_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_PROTEIN_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_LINKS_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_SVS_FILE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
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
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;
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
    public final List<VisGeneExon> Exons;
    public final List<Integer> Clusters;

    public final List<String> Chromosomes;
    public final Set<String> Genes;
    public final List<ChrBaseRegion> SpecificRegions;

    private final String mSampleDataDir;

    private static final String SAMPLE = "sample";
    private static final String VIS_FILE_DIRECTORY = "vis_file_dir";
    private static final String LOAD_COHORT_FILES = "load_cohort_files";
    private static final String CLUSTERS = "clusterId";
    private static final String CHROMOSOMES = "chromosome";
    private static final String GENE = "gene";

    public SampleData(final CommandLine cmd) throws Exception
    {
        final StringJoiner missingJoiner = new StringJoiner(", ");

        Sample = parameter(cmd, SAMPLE, missingJoiner);
        mSampleDataDir = cmd.getOptionValue(VIS_FILE_DIRECTORY);

        SpecificRegions = ChrBaseRegion.loadSpecificRegions(cmd);

        boolean useCohortFiles = cmd.hasOption(LOAD_COHORT_FILES);
        final String svDataFile = useCohortFiles ? COHORT_VIS_SVS_FILE : VisSvDataFile.generateFilename(mSampleDataDir, Sample);
        final String linksFile = useCohortFiles ? COHORT_VIS_LINKS_FILE : VisSegmentFile.generateFilename(mSampleDataDir, Sample);
        final String cnaFile = useCohortFiles ? COHORT_VIS_COPY_NUMBER_FILE : VisCopyNumberFile.generateFilename(mSampleDataDir, Sample);
        final String geneExonFile = useCohortFiles ? COHORT_VIS_GENE_EXONS_FILE : VisGeneExon.generateFilename(mSampleDataDir, Sample);
        final String proteinFile = useCohortFiles ? COHORT_VIS_PROTEIN_FILE : VisProteinDomainFile.generateFilename(mSampleDataDir, Sample);
        final String fusionFile = useCohortFiles ? COHORT_VIS_FUSIONS_FILE : VisFusionFile.generateFilename(mSampleDataDir, Sample);

        SvData = VisLinks.readSvData(svDataFile, Sample, SpecificRegions);

        Fusions = loadFusions(fusionFile).stream().filter(x -> x.sampleId().equals(Sample)).collect(toList());
        Exons = VisExons.readExons(geneExonFile).stream().filter(x -> x.SampleId.equals(Sample)).collect(toList());
        Segments = VisSegments.readTracks(linksFile).stream().filter(x -> x.sampleId().equals(Sample)).collect(toList());

        CopyNumberAlterations = VisCopyNumbers.read(cnaFile)
                .stream().filter(x -> x.sampleId().equals(Sample)).collect(toList());

        ProteinDomains = VisProteinDomains.readProteinDomains(proteinFile, Fusions)
                .stream()
                .filter(x -> x.sampleId().equals(Sample))
                .collect(toList());

        Clusters = parseClusters(cmd);
        Chromosomes = parseChromosomes(cmd);
        Genes = Sets.newHashSet();

        if(Segments.isEmpty() || SvData.isEmpty() || CopyNumberAlterations.isEmpty())
        {
            VIS_LOGGER.warn("sample({}) empty segments, SVs or copy-number files", Sample);
            return;
        }

        if(cmd.hasOption(GENE))
        {
            String geneStr = cmd.getOptionValue(GENE);
            String delim = geneStr.contains(";") ? ";" : ",";
            Arrays.stream(geneStr.split(delim)).forEach(x -> Genes.add(x));
        }

        boolean loadSvData = cmd.hasOption(PLOT_CLUSTER_GENES) || cmd.hasOption(RESTRICT_CLUSTERS_BY_GENE) || !SpecificRegions.isEmpty();

        if(loadSvData)
        {
            final String svAnnotationsFile = LinxSvAnnotation.generateFilename(mSampleDataDir, Sample);
            List<LinxSvAnnotation> svAnnotations = LinxSvAnnotation.read(svAnnotationsFile);

            if(!SpecificRegions.isEmpty())
            {
                svAnnotations = svAnnotations.stream().filter(x -> SvData.stream().anyMatch(y -> y.svId() == x.svId())).collect(toList());
            }

            if(cmd.hasOption(PLOT_CLUSTER_GENES) && !Clusters.isEmpty())
            {
                for(LinxSvAnnotation svAnnotation : svAnnotations)
                {
                    if(Clusters.contains(svAnnotation.clusterId()))
                    {
                        if(!svAnnotation.geneStart().isEmpty())
                            Arrays.stream(svAnnotation.geneStart().split(ITEM_DELIM)).forEach(x -> Genes.add(x));

                        if(!svAnnotation.geneEnd().isEmpty())
                            Arrays.stream(svAnnotation.geneEnd().split(ITEM_DELIM)).forEach(x -> Genes.add(x));
                    }
                }
            }
            else if(!Genes.isEmpty() && Clusters.isEmpty() && cmd.hasOption(RESTRICT_CLUSTERS_BY_GENE))
            {
                for(LinxSvAnnotation svAnnotation : svAnnotations)
                {
                    if(Genes.stream().anyMatch(x -> x.equals(svAnnotation.geneStart()) || x.equals(svAnnotation.geneEnd())))
                    {
                        if(!Clusters.contains(svAnnotation.clusterId()))
                            Clusters.add(svAnnotation.clusterId());
                    }
                }

                List<LinxDriver> linxDrivers = LinxDriver.read(LinxDriver.generateFilename(mSampleDataDir, Sample));

                for(LinxDriver linxDriver : linxDrivers)
                {
                    if(Genes.stream().anyMatch(x -> x.equals(linxDriver.gene())))
                    {
                        if(!Clusters.contains(linxDriver.clusterId()))
                            Clusters.add(linxDriver.clusterId());
                    }
                }
            }

            if(!SpecificRegions.isEmpty() && Clusters.isEmpty())
            {
                // limit to those clusters covered by an SV in the specific regions
                for(LinxSvAnnotation svAnnotation : svAnnotations)
                {
                    if(!Clusters.contains(svAnnotation.clusterId()))
                        Clusters.add(svAnnotation.clusterId());
                }
            }
        }

        Exons.addAll(additionalExons(Genes, cmd, Exons, Clusters));
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

    private static List<VisGeneExon> additionalExons(
            final Set<String> geneList, final CommandLine cmd, final List<VisGeneExon> currentExons, final List<Integer> clusterIds)
    {
        final List<VisGeneExon> exonList = Lists.newArrayList();

        if(geneList.isEmpty())
            return exonList;

        if(!cmd.hasOption(ENSEMBL_DATA_DIR))
            return exonList;

        final String sampleId = cmd.getOptionValue(SAMPLE);
        final List<Integer> allClusterIds = clusterIds.isEmpty() ? Lists.newArrayList(0) : clusterIds;

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, V37.toString()));

        EnsemblDataCache geneTransCache = new EnsemblDataCache(cmd, refGenomeVersion);
        geneTransCache.setRequiredData(true, false, false, true);
        geneTransCache.load(false);

        for(final String geneName : geneList)
        {
            if(currentExons.stream().anyMatch(x -> x.Gene.equals(geneName)))
                continue;

            VIS_LOGGER.info("loading exon data for additional gene({})", geneName);

            GeneData geneData = geneTransCache.getGeneDataByName(geneName);
            TranscriptData transcriptData = geneData != null ? geneTransCache.getCanonicalTranscriptData(geneData.GeneId) : null;

            if(transcriptData == null)
            {
                VIS_LOGGER.warn("data not found for specified gene({})", geneName);
                continue;
            }

            for(Integer clusterId : allClusterIds)
            {
                exonList.addAll(VisExons.extractExonList(sampleId, clusterId, geneData, transcriptData));
            }
        }

        return exonList;
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(SAMPLE, true, "Sample name");
        options.addOption(VIS_FILE_DIRECTORY, true, "Path to all Linx vis files, used instead of specifying them individually");
        options.addOption(LOAD_COHORT_FILES, false, "Load Linx cohort rather than per-sample vis files");
        options.addOption(GENE, true, "Show canonical transcript for genes (separated by ','");
        options.addOption(RESTRICT_CLUSTERS_BY_GENE, false, "Only plot clusters with breakends in configured 'gene' list");
        addEnsemblDir(options);
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(CLUSTERS, true, "Only generate image for specified comma separated clusters");
        options.addOption(CHROMOSOMES, true, "Only generate image for specified comma separated chromosomes");
        options.addOption(SPECIFIC_REGIONS, true, SPECIFIC_REGIONS_DESC);
    }
}
