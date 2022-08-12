package com.hartwig.hmftools.geneutils.drivers;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
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
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

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
    private final List<String> mPanelGeneOverrides;

    // config
    private static final String DRIVER_GENE_PANEL_TSV = "driver_gene_panel";
    private static final String RESOURCE_REPO_DIR = "resource_repo_dir";

    private static final String GENE_PANEL_DIR = "gene_panel";
    private static final String SAGE_DIR = "sage";
    private static final String ENSEMBL_DIR = "ensembl_data_cache";
    private static final String PANEL_GENE_OVERRIDES = "panel_gene_overrides";

    public GenerateDriverGeneFiles(final CommandLine cmd)
    {
        GU_LOGGER.info("starting driver gene panel generation");

        mDriverGenePanelFile = cmd.getOptionValue(DRIVER_GENE_PANEL_TSV);
        mResourceRepoDir = checkAddDirSeparator(cmd.getOptionValue(RESOURCE_REPO_DIR));
        mOutputDir = parseOutputDir(cmd);

        mPanelGeneOverrides = Lists.newArrayList();

        if(cmd.hasOption(PANEL_GENE_OVERRIDES))
        {
            Arrays.stream(cmd.getOptionValue(PANEL_GENE_OVERRIDES).split(",", -1)).forEach(x -> mPanelGeneOverrides.add(x));
        }
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

        writeGermlineBlacklist(refGenomeVersion, sageDir);

        writeGermlineHotspots(refGenomeVersion, driverGenes, sageDir);
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

        List<NamedBed> panelRegionsWithUtr = Lists.newArrayList();
        List<NamedBed> panelRegionsWithoutUtr = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            List<GeneData> geneDataList = chrGeneDataMap.get(refGenomeVersion.versionedChromosome(chromosome.toString()));

            for(GeneData geneData : geneDataList)
            {
                DriverGene driverGene = driverGenes.stream().filter(x -> x.gene().equals(geneData.GeneName)).findFirst().orElse(null);

                if(driverGene == null)
                    continue;

                if(!driverGene.reportSomatic() && !driverGene.reportGermline() && !mPanelGeneOverrides.contains(geneData.GeneName))
                    continue;

                TranscriptData transData = ensemblDataCache.getTranscriptData(geneData.GeneId, "");

                List<GenomeRegion> regionsWithUtr = getTranscriptRegions(geneData, transData, true);
                List<GenomeRegion> regionsWithoutUtr = getTranscriptRegions(geneData, transData, false);

                // merge any overlap with the previous gene region
                for(int i = 0; i <= 1; ++i)
                {
                    List<GenomeRegion> newRegions = i == 0 ? regionsWithUtr : regionsWithoutUtr;

                    if(newRegions.isEmpty())
                        continue;

                    List<NamedBed> allRegions = i == 0 ? panelRegionsWithUtr : panelRegionsWithoutUtr;

                    if(allRegions.isEmpty())
                        continue;

                    NamedBed lastRegion = allRegions.get(allRegions.size() - 1);
                    if(!lastRegion.chromosome().equals(chromosome.toString()))
                        continue;

                    int newLastRegionEnd = 0;
                    int regionsRemoved = 0;
                    while(!newRegions.isEmpty())
                    {
                        GenomeRegion newRegion = newRegions.get(0);
                        if(newRegion.start() > lastRegion.end() + 1)
                            break;

                        // otherwise remove the new region and merge
                        newLastRegionEnd = max(newRegion.end(), lastRegion.end());
                        newRegions.remove(0);
                        ++regionsRemoved;
                    }

                    if(newLastRegionEnd > 0)
                    {
                        NamedBed newLastRegion = ImmutableNamedBed.builder().from(lastRegion).end(newLastRegionEnd).build();
                        allRegions.set(allRegions.size() - 1, newLastRegion);

                        GU_LOGGER.debug("gene({}) removed {} regions from overlap with previous region(gene={} range={}->{})",
                                geneData.GeneName, regionsRemoved, lastRegion.name(), lastRegion.start(), lastRegion.end());
                    }
                }

                regionsWithUtr.stream()
                        .map(x -> ImmutableNamedBed.builder().from(x).name(geneData.GeneName).build())
                        .forEach(x -> panelRegionsWithUtr.add(x));

                regionsWithoutUtr.stream()
                        .map(x -> ImmutableNamedBed.builder().from(x).name(geneData.GeneName).build())
                        .forEach(x -> panelRegionsWithoutUtr.add(x));
            }
        }

        String codingWithUtr = formVersionFile(sageDir, "ActionableCodingPanel.bed.gz", refGenomeVersion);
        GU_LOGGER.info("writing {} panel coding regions file({})", refGenomeVersion, codingWithUtr);

        String coverageWithoutUtr = formVersionFile(sageDir, "CoverageCodingPanel.bed.gz", refGenomeVersion);
        GU_LOGGER.info("writing {} panel coverage regions file({})", refGenomeVersion, coverageWithoutUtr);

        try
        {
            NamedBedFile.writeBedFile(codingWithUtr, panelRegionsWithUtr);
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

            if(positionsOverlap(startPosition, endPosition, exonStart, exonEnd))
            {
                regions.add(GenomeRegions.create(geneData.Chromosome, max(startPosition, exonStart), min(endPosition, exonEnd)));
            }
        }

        return regions;
    }

    private void writeGermlineBlacklist(final RefGenomeVersion refGenomeVersion, final String sageDir)
    {
        String germlineBlacklistFile = formVersionFile(sageDir, "KnownBlacklist.germline.vcf.gz", refGenomeVersion);

        GU_LOGGER.info("writing {} germline blacklist file at {}", refGenomeVersion, germlineBlacklistFile);

        try
        {
            List<VariantContext> germlineBlackList = refGenomeVersion == RefGenomeVersion.V37 ?
                    GermlineResources.blacklist37() : GermlineResources.blacklist38();

            GermlineBlacklistVCF.write(germlineBlacklistFile, germlineBlackList);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write germline blacklist file: {}", e.toString());
        }
    }

    private void writeGermlineHotspots(final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes, final String sageDir)
    {
        String germlineHotspotFile = formVersionFile(sageDir, "KnownHotspots.germline.vcf.gz", refGenomeVersion);

        String sageRefDir = mResourceRepoDir + SAGE_DIR + File.separator + refGenomeVersion.identifier();
        String clinvarFile = formVersionFile(sageRefDir, "clinvar.vcf.gz", refGenomeVersion);

        GU_LOGGER.info("located clinvar file for {} at {}", refGenomeVersion, clinvarFile);

        List<String> germlineHotspotGenes = driverGenes.stream()
                .filter(x -> x.reportGermlineHotspot() != DriverGeneGermlineReporting.NONE)
                .map(x -> x.gene()).collect(Collectors.toList());

        GU_LOGGER.info("writing {} germline hotspots to {}", germlineHotspotGenes.size(), germlineHotspotFile);

        try
        {
            GermlineHotspotVCF.write(refGenomeVersion, clinvarFile, germlineHotspotFile, germlineHotspotGenes);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write germline hotspots file: {}", e.toString());
        }
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

    public static void main(String[] args) throws IOException, ParseException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        setLogLevel(cmd);

        GenerateDriverGeneFiles generator = new GenerateDriverGeneFiles(cmd);
        generator.run();
    }

    private static Options createOptions()
    {
        Options options = new Options();

        options.addOption(DRIVER_GENE_PANEL_TSV, true, "File containing the driver gene panel for 37");
        options.addOption(RESOURCE_REPO_DIR, true, "The directory holding the public hmf resources repo");
        options.addOption(PANEL_GENE_OVERRIDES, true, "List of comma-separated genes to include in panel");
        addOutputDir(options);
        addLoggingOptions(options);

        return options;
    }
}
