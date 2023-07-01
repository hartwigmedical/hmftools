package com.hartwig.hmftools.geneutils.drivers;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.RESOURCE_REPO_DIR;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.RESOURCE_REPO_DIR_DESC;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.getEnsemblDirectory;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.logVersion;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.variant.variantcontext.VariantContext;

public class GenerateDriverGeneFiles
{
    private final String mDriverGenePanelFile;
    private final String mResourceRepoDir;
    private final String mOutputDir;
    private final List<String> mPanelGeneOverrides;

    // config
    private static final String DRIVER_GENE_PANEL_TSV = "driver_gene_panel";

    private static final String GENE_PANEL_DIR = "gene_panel";
    private static final String SAGE_DIR = "sage";
    private static final String PANEL_GENE_OVERRIDES = "panel_gene_overrides";

    public GenerateDriverGeneFiles(final ConfigBuilder configBuilder)
    {
        GU_LOGGER.info("starting driver gene panel generation");

        mDriverGenePanelFile = configBuilder.getValue(DRIVER_GENE_PANEL_TSV);
        mResourceRepoDir = checkAddDirSeparator(configBuilder.getValue(RESOURCE_REPO_DIR));
        mOutputDir = parseOutputDir(configBuilder);

        mPanelGeneOverrides = Lists.newArrayList();

        if(configBuilder.hasValue(PANEL_GENE_OVERRIDES))
        {
            Arrays.stream(configBuilder.getValue(PANEL_GENE_OVERRIDES).split(",", -1)).forEach(x -> mPanelGeneOverrides.add(x));
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
        return format("%s/%s", outputDir, refGenomeVersion.addVersionToFilePath(filename));
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
        String ensemblDir = getEnsemblDirectory(refGenomeVersion, mResourceRepoDir);

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);

        Set<String> actionablePanel = driverGenes.stream()
                .filter(x -> x.reportSomatic() || x.reportGermline() || x.reportPGX() || mPanelGeneOverrides.contains(x.gene()))
                .map(x -> x.gene())
                .collect(Collectors.toSet());

        String codingWithUtr = formVersionFile(sageDir, "ActionableCodingPanel.bed.gz", refGenomeVersion);

        GU_LOGGER.info("writing {} panel coding regions file({}) for {} genes",
                refGenomeVersion, codingWithUtr, actionablePanel.size());

        writeGenePanelRegions(refGenomeVersion, ensemblDataCache, actionablePanel, true, codingWithUtr);

        Set<String> coveragePanel = driverGenes.stream().map(x -> x.gene()).collect(Collectors.toSet());

        String coverageWithoutUtr = formVersionFile(sageDir, "CoverageCodingPanel.bed.gz", refGenomeVersion);

        GU_LOGGER.info("writing {} panel coverage regions file({}) for {} genes",
                refGenomeVersion, coverageWithoutUtr, coveragePanel.size());

        writeGenePanelRegions(refGenomeVersion, ensemblDataCache, coveragePanel, false, coverageWithoutUtr);
    }

    private void writeGenePanelRegions(
            final RefGenomeVersion refGenomeVersion, final EnsemblDataCache ensemblDataCache,
            final Set<String> geneSet, boolean includeUTR, final String outputFile)
    {
        final Map<String,List<GeneData>> chrGeneDataMap = ensemblDataCache.getChrGeneDataMap();

        List<CodingRegion> panelRegions = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = refGenomeVersion.versionedChromosome(chromosome.toString());
            List<GeneData> geneDataList = chrGeneDataMap.get(chromosomeStr);

            for(GeneData geneData : geneDataList)
            {
                if(!geneSet.contains(geneData.GeneName))
                    continue;

                TranscriptData transData = ensemblDataCache.getTranscriptData(geneData.GeneId, "");

                List<CodingRegion> transcriptRegions = getTranscriptRegions(geneData, transData, includeUTR);

                // merge any overlap with the previous gene region
                CodingRegion lastRegion = !panelRegions.isEmpty() ? panelRegions.get(panelRegions.size() - 1) : null;

                if(!transcriptRegions.isEmpty() && lastRegion != null && lastRegion.Chromosome.equals(chromosomeStr))
                {
                    int newLastRegionEnd = 0;
                    int regionsRemoved = 0;
                    while(!transcriptRegions.isEmpty())
                    {
                        CodingRegion newRegion = transcriptRegions.get(0);
                        if(newRegion.start() > lastRegion.end() + 1)
                            break;

                        // otherwise remove the new region and merge
                        newLastRegionEnd = max(newRegion.end(), lastRegion.end());
                        transcriptRegions.remove(0);
                        ++regionsRemoved;
                    }

                    if(newLastRegionEnd > 0)
                    {
                        CodingRegion newLastRegion = new CodingRegion(
                                lastRegion.Chromosome, lastRegion.start(), newLastRegionEnd, lastRegion.GeneName, lastRegion.ExonRank);

                        panelRegions.set(panelRegions.size() - 1, newLastRegion);

                        GU_LOGGER.debug("gene({}) removed {} regions from overlap with previous region(gene={} range={}->{})",
                                geneData.GeneName, regionsRemoved, lastRegion.GeneName, lastRegion.start(), lastRegion.end());
                    }
                }

                panelRegions.addAll(transcriptRegions);
            }
        }

        try
        {
            List<NamedBed> bedRegions = panelRegions.stream().map(x -> x.asBed()).collect(Collectors.toList());

            NamedBedFile.writeBedFile(outputFile, bedRegions);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write gene panel file: {}", e.toString());
        }
    }

    private class CodingRegion extends ChrBaseRegion
    {
        public final String GeneName;
        public final int ExonRank;

        public CodingRegion(final String chromosome, final int posStart, final int posEnd, final String geneName, final int exonRank)
        {
            super(chromosome, posStart, posEnd);
            GeneName = geneName;
            ExonRank = exonRank;
        }

        public NamedBed asBed()
        {
            return ImmutableNamedBed.builder()
                    .chromosome(Chromosome)
                    .start(start())
                    .end(end())
                    .name(format("%s_%d", GeneName, ExonRank))
                    .build();
        }
    }

    private static final int SPLICE_SIZE = 10;

    private List<CodingRegion> getTranscriptRegions(final GeneData geneData, final TranscriptData transData, boolean includeUTR)
    {
        int startPosition = includeUTR || transData.nonCoding() ? transData.TransStart : transData.CodingStart;
        int endPosition = includeUTR || transData.nonCoding() ? transData.TransEnd : transData.CodingEnd;

        final List<CodingRegion> regions = Lists.newArrayList();

        for(int i = 0; i < transData.exons().size(); i++)
        {
            ExonData exon = transData.exons().get(i);
            int exonStart = i == 0 ? exon.Start : exon.Start - SPLICE_SIZE;
            int exonEnd = i == transData.exons().size() - 1 ? exon.End : exon.End + SPLICE_SIZE;

            if(positionsOverlap(startPosition, endPosition, exonStart, exonEnd))
            {
                regions.add(new CodingRegion(
                        geneData.Chromosome, max(startPosition, exonStart), min(endPosition, exonEnd), geneData.GeneName, exon.Rank));
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

    public static void main(String[] args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        addGenePanelOption(configBuilder, true);
        configBuilder.addPath(RESOURCE_REPO_DIR, true, RESOURCE_REPO_DIR_DESC);
        configBuilder.addConfigItem(PANEL_GENE_OVERRIDES, "List of comma-separated genes to include in panel");
        addOutputDir(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);
        logVersion();

        GenerateDriverGeneFiles generator = new GenerateDriverGeneFiles(configBuilder);
        generator.run();
    }
}
