package com.hartwig.hmftools.geneutils.drivers;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.bed.ImmutableNamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.genepanel.HmfExonPanelBed;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.variant.variantcontext.VariantContext;

public class GenerateDriverGeneFiles
{
    private final String mDriverGenePanelFile;
    private final String mResourceRepoDir;
    private final String mOutputDir;

    // config
    private static final String DRIVER_GENE_PANEL_TSV = "driver_gene_panel";
    private static final String RESOURCE_REPO_DIR = "resource_repo_dir";

    private static final String GENE_PANEL_DIR = "gene_panel";
    private static final String SAGE_DIR = "sage";
    private static final String ENSEMBL_DIR = "ensembl_data_cache";

    public GenerateDriverGeneFiles(final CommandLine cmd)
    {
        GU_LOGGER.info("starting driver gene panel generation");

        mDriverGenePanelFile = cmd.getOptionValue(DRIVER_GENE_PANEL_TSV);
        mResourceRepoDir = checkAddDirSeparator(cmd.getOptionValue(RESOURCE_REPO_DIR));
        mOutputDir = parseOutputDir(cmd);
    }

    public void run() throws IOException
    {
        GU_LOGGER.info("resource reference directory: {}", mResourceRepoDir);
        GU_LOGGER.info("output directory: {}", mOutputDir);

        createOutputDir(mOutputDir);
        createOutputDir(mOutputDir + GENE_PANEL_DIR + File.separator);
        createOutputDir(mOutputDir + SAGE_DIR + File.separator);

        List<DriverGene> driverGenes = DriverGeneFile.read(mDriverGenePanelFile);
        GU_LOGGER.info("loaded {} driver genes from {}", driverGenes.size(), mDriverGenePanelFile);

        process(RefGenomeVersion.V37, driverGenes);
        process(RefGenomeVersion.V38, driverGenes);

        GU_LOGGER.info("file generation complete");
    }

    private static String formVersionFile(
            final String outputDir, final String filename, final RefGenomeVersion refGenomeVersion)
    {
        return String.format("%s/%s", outputDir, refGenomeVersion.addVersionToFilePath(filename));
    }

    public void process(final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes) throws IOException
    {
        Collections.sort(driverGenes);
        writeDriverGeneFiles(refGenomeVersion, driverGenes);

        String sageDir = mOutputDir + SAGE_DIR + File.separator + refGenomeVersion.identifier();
        createOutputDir(sageDir + File.separator);

        writeGenePanelRegions(refGenomeVersion, driverGenes, sageDir);

        String germlineBlacklistFile = formVersionFile(sageDir, "KnownBlacklist.germline.vcf.gz", refGenomeVersion);

        GU_LOGGER.info("writing {} germline blacklist file at {}", refGenomeVersion, germlineBlacklistFile);
        List<VariantContext> germlineBlackList = refGenomeVersion == RefGenomeVersion.V37 ?
                GermlineResources.blacklist37() : GermlineResources.blacklist38();

        GermlineBlacklistVCF.process(germlineBlacklistFile, germlineBlackList);

        String germlineHotspotFile = formVersionFile(sageDir, "KnownHotspots.germline.vcf.gz", refGenomeVersion);

        GU_LOGGER.info("writing {} germline hotspots at {}", refGenomeVersion, germlineHotspotFile);
        GermlineHotspotVCF germlineHotspotVCF = new GermlineHotspotVCF(germlineHotspotGenes(driverGenes));

        /*
        String somaticCodingWithoutUtr = formVersionFile(sageDir, "ActionableCodingPanel.somatic.bed.gz", refGenomeVersion);
        String germlineCodingWithUtr = formVersionFile(sageDir, "ActionableCodingPanel.germline.bed.gz", refGenomeVersion);
        String germlineCodingWithoutUtr = formVersionFile(sageDir, "CoverageCodingPanel.germline.bed.gz", refGenomeVersion);
        String germlineSliceFile = formVersionFile(sageDir, "SlicePanel.germline.bed.gz", refGenomeVersion);

        Set<String> germlineGenes = germlineGenes(driverGenes);
        Set<String> somaticGenes = somaticGenes(driverGenes);
        Set<String> allGenes = Sets.newHashSet();
        allGenes.addAll(germlineGenes);
        allGenes.addAll(somaticGenes);

        // Write out actionable bed files
        List<HmfTranscriptRegion> transcripts = refGenomeVersion == RefGenomeVersion.V37 ?
                HmfGenePanelSupplier.allGeneList37() : HmfGenePanelSupplier.allGeneList38();

        GU_LOGGER.info("writing {} {} bed file at {}", refGenomeVersion, "somatic", somaticCodingWithoutUtr);
        createNamedBedFiles(false, somaticCodingWithoutUtr, somaticGenes, transcripts);

        GU_LOGGER.info("writing {} {} coverage bed file at {}", refGenomeVersion, "germline", germlineCodingWithoutUtr);
        createNamedBedFiles(false, germlineCodingWithoutUtr, allGenes, transcripts);

        GU_LOGGER.info("writing {} {} bed file at {}", refGenomeVersion, "germline", germlineCodingWithUtr);
        List<GenomeRegion> germlinePanel = createUnnamedBedFiles(true, germlineCodingWithUtr, allGenes, transcripts);

        String sageRefDir = mResourceRepoDir + SAGE_DIR + File.separator + refGenomeVersion.identifier();

        String clinvarFile = formVersionFile(sageRefDir, "clinvar.vcf.gz", refGenomeVersion);

        GU_LOGGER.info("located clinvar file for {} at {}", refGenomeVersion, clinvarFile);

        List<GenomeRegion> germlineHotspots = germlineHotspotVCF.process(clinvarFile, germlineHotspotFile);

        String qualityBedFile = getResourceURL(refGenomeVersion.addVersionToFilePath("/drivers/QualityRecalibration.bed"));
        Collection<GenomeRegion> qualityRecalibrationRegions = BEDFileLoader.fromBedFile(qualityBedFile).values();

        GU_LOGGER.info("writing {} germline slice file at {}", refGenomeVersion, germlineSliceFile);
        createSliceFile(germlineSliceFile, germlinePanel, germlineHotspots, qualityRecalibrationRegions);
        */
    }

    private void writeDriverGeneFiles(final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes)
    {
        try
        {
            String genePanelDir = mOutputDir + GENE_PANEL_DIR + File.separator + refGenomeVersion.identifier() + File.separator;
            createOutputDir(genePanelDir);
            String driverGeneFile = refGenomeVersion.addVersionToFilePath(genePanelDir + "DriverGenePanel.tsv");
            DriverGeneFile.write(driverGeneFile, driverGenes);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write driver gene panel files: {}", e.toString());
        }
    }

    private void writeGenePanelRegions(final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes, final String sageDir)
    {
        String ensemblDir = mResourceRepoDir + ENSEMBL_DIR + File.separator + refGenomeVersion.identifier();

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);

        final Map<String,List<GeneData>> chrGeneDataMap = ensemblDataCache.getChrGeneDataMap();

        List<GenomeRegion> panelRegionsWithUtr = Lists.newArrayList();
        List<NamedBed> panelRegionsWithoutUtr = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            List<GeneData> geneDataList = chrGeneDataMap.get(refGenomeVersion.versionedChromosome(chromosome.toString()));

            for(GeneData geneData : geneDataList)
            {
                if(!driverGenes.stream().anyMatch(x -> x.gene().equals(geneData.GeneName)))
                    continue;

                TranscriptData transData = ensemblDataCache.getTranscriptData(geneData.GeneId, "");

                panelRegionsWithUtr.addAll(getTranscriptRegions(geneData, transData, true));

                List<GenomeRegion> regions = getTranscriptRegions(geneData, transData, false);

                regions.stream()
                        .map(x -> ImmutableNamedBed.builder().from(x).name(geneData.GeneName).build())
                        .forEach(x -> panelRegionsWithoutUtr.add(x));
            }
        }

        String codingWithUtr = formVersionFile(sageDir, "ActionableCodingPanel.bed.gz", refGenomeVersion);
        String coverageWithoutUtr = formVersionFile(sageDir, "CoverageCodingPanel.bed.gz", refGenomeVersion);

        GU_LOGGER.info("writing {} panel coverage regions file({})", refGenomeVersion, coverageWithoutUtr);
        GU_LOGGER.info("writing {} panel coding regions file({})", refGenomeVersion, codingWithUtr);

        try
        {
            NamedBedFile.writeUnnamedBedFile(codingWithUtr, panelRegionsWithUtr);
            NamedBedFile.writeBedFile(coverageWithoutUtr, panelRegionsWithoutUtr);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write gene panel files: {}", e.toString());
        }
    }

    private static final int SPLICE_SIZE = 10;

    private List<GenomeRegion> getTranscriptRegions(final GeneData geneData, final TranscriptData transData, boolean includeUTR)
    {
        int startPosition = includeUTR || transData.nonCoding() ? transData.TransStart : transData.CodingStart;
        int endPosition = includeUTR || transData.nonCoding() ? transData.TransEnd : transData.CodingEnd;

        final List<GenomeRegion> regions = Lists.newArrayList();

        for(int i = 0; i < transData.exons().size(); i++)
        {
            ExonData exon = transData.exons().get(i);
            int exonStart = i == 0 ? exon.Start : exon.Start - SPLICE_SIZE;
            int exonEnd = i == transData.exons().size() - 1 ? exon.End : exon.End + SPLICE_SIZE;

            if(startPosition < exonEnd && endPosition > exonStart)
            {
                regions.add(GenomeRegions.create(geneData.Chromosome, max(startPosition, exonStart), min(endPosition, exonEnd)));
            }
        }

        return regions;
    }

    private boolean createOutputDir(final String outputDir)
    {
        final File dir = new File(outputDir);
        if(!dir.exists() && !dir.mkdirs())
        {
            GU_LOGGER.error("unable to write directory " + outputDir);
            return false;
        }

        return true;
    }

    private static void createNamedBedFiles(
            boolean includeUTR, final String file, final Set<String> genes, final List<HmfTranscriptRegion> transcripts)
    {
        try
        {
            List<NamedBed> somaticBed = HmfExonPanelBed.createNamedCodingRegions(includeUTR, genes, transcripts);
            NamedBedFile.writeBedFile(file, somaticBed);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write name bed file({}): {}", file, e.toString());
        }
    }

    private static void createSliceFile(final String file, final Collection<GenomeRegion>... allRegions) throws IOException
    {
        GenomeRegionsBuilder builder = new GenomeRegionsBuilder();
        for(Collection<GenomeRegion> regions : allRegions)
        {
            for(GenomeRegion region : regions)
            {
                builder.addRegion(region.chromosome(), region.start() - 100, region.end() + 100);
            }
        }

        List<GenomeRegion> combined = builder.build();
        NamedBedFile.writeUnnamedBedFile(file, combined);
    }


    private static List<GenomeRegion> createUnnamedBedFiles(
            boolean includeUTR, final String file, final Set<String> genes,
            final List<HmfTranscriptRegion> transcripts) throws IOException
    {
        List<GenomeRegion> somaticBed = HmfExonPanelBed.createUnnamedCodingRegions(includeUTR, genes, transcripts);
        NamedBedFile.writeUnnamedBedFile(file, somaticBed);
        return somaticBed;
    }

    private static Set<String> somaticGenes(final List<DriverGene> genePanel)
    {
        Set<String> somaticReportableGenes = Sets.newHashSet();
        for(DriverGene driverGene : genePanel)
        {
            if(driverGene.reportSomatic())
            {
                somaticReportableGenes.add(driverGene.gene());
            }
        }

        return somaticReportableGenes;
    }

    private static Set<String> germlineGenes(final List<DriverGene> genePanel)
    {
        Set<String> germlineReportableGenes = Sets.newHashSet();
        for(DriverGene driverGene : genePanel)
        {
            if(driverGene.reportGermline())
            {
                germlineReportableGenes.add(driverGene.gene());
            }
        }

        return germlineReportableGenes;
    }

    private static Set<String> germlineHotspotGenes(final List<DriverGene> genePanel)
    {
        Set<String> germlineHotspotGenes = Sets.newHashSet();
        for(DriverGene driverGene : genePanel)
        {
            if(driverGene.reportGermlineHotspot() != DriverGeneGermlineReporting.NONE)
            {
                germlineHotspotGenes.add(driverGene.gene());
            }
        }

        return germlineHotspotGenes;
    }

    private static String getResourceURL(final String location)
    {
        return GenerateDriverGeneFiles.class.getResource(location).toString();
    }

    public static void main(String[] args) throws IOException, ParseException
    {
        // See also https://github.com/hartwigmedical/scratchpad/tree/master/genePanel
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        GenerateDriverGeneFiles generator = new GenerateDriverGeneFiles(cmd);
        generator.run();
    }

    private static Options createOptions()
    {
        Options options = new Options();

        options.addOption(DRIVER_GENE_PANEL_TSV, true, "File containing the driver gene panel for 37");
        options.addOption(RESOURCE_REPO_DIR, true, "The directory holding the public hmf resources repo");
        addOutputDir(options);
        addLoggingOptions(options);

        return options;
    }
}
