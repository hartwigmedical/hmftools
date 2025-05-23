package com.hartwig.hmftools.geneutils.drivers;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneRegions.getTranscriptRegions;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.RESOURCE_REPO_DIR;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.RESOURCE_REPO_DIR_DESC;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.createOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.getEnsemblDirectory;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneRegion;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.variant.variantcontext.VariantContext;

public class GenerateDriverGeneFiles
{
    private final String mDriverGenePanelFile;
    private final String mGeneIdFile;
    private final String mResourceRepoDir;
    private final String mOutputDir;
    private final List<String> mPanelGeneOverrides;
    private final List<RefGenomeVersion> mRefGenomeVersions;

    private static final String SAGE_RESOURCE_DIR = "sage";

    // config
    private static final String PANEL_GENE_OVERRIDES = "panel_gene_overrides";

    public GenerateDriverGeneFiles(final ConfigBuilder configBuilder)
    {
        GU_LOGGER.info("starting driver gene panel generation");

        mDriverGenePanelFile = configBuilder.getValue(DRIVER_GENE_PANEL);
        mGeneIdFile = configBuilder.getValue(GENE_ID_FILE);
        mResourceRepoDir = checkAddDirSeparator(configBuilder.getValue(RESOURCE_REPO_DIR));
        mOutputDir = parseOutputDir(configBuilder);

        mRefGenomeVersions = Lists.newArrayList();
        mPanelGeneOverrides = Lists.newArrayList();

        if(configBuilder.hasValue(PANEL_GENE_OVERRIDES))
        {
            Arrays.stream(configBuilder.getValue(PANEL_GENE_OVERRIDES).split(",", -1)).forEach(x -> mPanelGeneOverrides.add(x));
        }

        if(configBuilder.hasValue(REF_GENOME_VERSION))
        {
            mRefGenomeVersions.add(RefGenomeVersion.from(configBuilder));
        }
        else
        {
            mRefGenomeVersions.add(V37);
            mRefGenomeVersions.add(V38);
        }
    }

    public void run() throws IOException
    {
        GU_LOGGER.info("resource reference directory: {}", mResourceRepoDir);
        GU_LOGGER.info("output directory: {}", mOutputDir);

        createOutputDir(mOutputDir);

        List<String> actionableGenes = Lists.newArrayList(mPanelGeneOverrides);
        List<String> coverageGenes = Lists.newArrayList(mPanelGeneOverrides);

        List<DriverGene> driverGenes = Lists.newArrayList();

        if(mDriverGenePanelFile != null)
        {
            driverGenes.addAll(DriverGeneFile.read(mDriverGenePanelFile));

            Collections.sort(driverGenes);

            GU_LOGGER.info("loaded {} driver genes from {}", driverGenes.size(), mDriverGenePanelFile);

            driverGenes.stream()
                    .filter(x -> !mPanelGeneOverrides.contains(x.gene()))
                    .filter(x -> x.reportSomatic() || x.reportGermline() || x.reportPGX())
                    .forEach(x -> actionableGenes.add(x.gene()));

            driverGenes.forEach(x -> coverageGenes.add(x.gene()));
        }
        else if(mGeneIdFile != null)
        {
            List<String> geneNames = loadDelimitedIdFile(mGeneIdFile, FLD_GENE_NAME, CSV_DELIM);

            Collections.sort(geneNames);

            actionableGenes.addAll(geneNames);
            coverageGenes.addAll(geneNames);
        }

        for(RefGenomeVersion refGenomeVersion : mRefGenomeVersions)
        {
            process(refGenomeVersion, driverGenes, actionableGenes, coverageGenes);
        }

        GU_LOGGER.info("file generation complete");
    }

    private static String formVersionFile(
            final String outputDir, final String filename, final RefGenomeVersion refGenomeVersion)
    {
        return format("%s/%s", outputDir, refGenomeVersion.addVersionToFilePath(filename));
    }

    public void process(
            final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes,
            final List<String> actionableGenes, final List<String> coverageGenes)
    {
        if(!driverGenes.isEmpty())
        {
            writeDriverGeneFiles(refGenomeVersion, driverGenes);

            writeGermlineBlacklist(refGenomeVersion);
            writeGermlineHotspots(refGenomeVersion, driverGenes);
        }

        writeGenePanelRegions(refGenomeVersion, actionableGenes, coverageGenes);
    }

    private void writeDriverGeneFiles(final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes)
    {
        try
        {
            String driverGeneFile = refGenomeVersion.addVersionToFilePath(mOutputDir + "DriverGenePanel.tsv");
            DriverGeneFile.write(driverGeneFile, driverGenes);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write driver gene panel files: {}", e.toString());
        }
    }

    private void writeGenePanelRegions(
            final RefGenomeVersion refGenomeVersion, final List<String> actionableGenes, final List<String> coverageGenes)
    {
        String ensemblDir = getEnsemblDirectory(refGenomeVersion, mResourceRepoDir);

        GU_LOGGER.debug("loading Ensembl data cache");

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);

        String codingWithUtr = formVersionFile(mOutputDir, "ActionableCodingPanel.bed.gz", refGenomeVersion);

        GU_LOGGER.info("writing {} panel coding regions file({}) for {} genes",
                refGenomeVersion, codingWithUtr, actionableGenes.size());

        writeGenePanelRegions(refGenomeVersion, ensemblDataCache, actionableGenes, true, codingWithUtr);
    }

    private void writeGenePanelRegions(
            final RefGenomeVersion refGenomeVersion, final EnsemblDataCache ensemblDataCache,
            final List<String> geneSet, boolean includeUTR, final String outputFile)
    {
        final Map<String,List<GeneData>> chrGeneDataMap = ensemblDataCache.getChrGeneDataMap();

        List<GeneRegion> panelRegions = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = refGenomeVersion.versionedChromosome(chromosome.toString());
            List<GeneData> geneDataList = chrGeneDataMap.get(chromosomeStr);

            List<GeneRegion> chrPanelRegions = Lists.newArrayList();

            for(GeneData geneData : geneDataList)
            {
                if(!geneSet.contains(geneData.GeneName))
                    continue;

                TranscriptData transData = ensemblDataCache.getTranscriptData(geneData.GeneId, "");

                List<GeneRegion> transcriptRegions = getTranscriptRegions(geneData, transData, includeUTR, true);

                chrPanelRegions.addAll(transcriptRegions);
            }

            // sort and merge any overlaps
            Collections.sort(chrPanelRegions);

            int regionsRemoved = 0;

            int index = 0;

            // check for overlaps with the previous region
            while(index < chrPanelRegions.size() - 1)
            {
                GeneRegion region = chrPanelRegions.get(index);

                int nextIndex = index + 1;
                while(nextIndex < chrPanelRegions.size())
                {
                    GeneRegion nextRegion = chrPanelRegions.get(nextIndex);

                    if(region.end() >= nextRegion.start())
                    {
                        GU_LOGGER.trace("gene({}) merged region({}) with next({})", region.GeneName, region, nextRegion);

                        if(nextRegion.end() > region.end())
                        {
                            region.setEnd(nextRegion.end());
                        }
                        ++regionsRemoved;
                        chrPanelRegions.remove(nextIndex);
                    }
                    else
                    {
                        break;
                    }
                }

                ++index;
            }

            if(regionsRemoved > 0)
            {
                GU_LOGGER.debug("chr({}) merged {} regions from overlaps", chromosomeStr, regionsRemoved);
            }

            panelRegions.addAll(chrPanelRegions);
        }

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            for(GeneRegion geneRegion : panelRegions)
            {
                String geneExonRank = format("%s_%d", geneRegion.GeneName, geneRegion.ExonRank);

                writer.write(format("%s\t%d\t%d\t%s",
                        geneRegion.Chromosome, geneRegion.start() - 1, geneRegion.end(), geneExonRank));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write gene panel file: {}", e.toString());
        }
    }

    private void writeGermlineBlacklist(final RefGenomeVersion refGenomeVersion)
    {
        String germlineBlacklistFile = formVersionFile(mOutputDir, "KnownBlacklist.germline.vcf.gz", refGenomeVersion);

        GU_LOGGER.info("writing {} germline blacklist file at {}", refGenomeVersion, germlineBlacklistFile);

        try
        {
            List<VariantContext> germlineBlackList = refGenomeVersion == V37 ?
                    GermlineResources.blacklist37() : GermlineResources.blacklist38();

            GermlineBlacklistVCF.write(germlineBlacklistFile, germlineBlackList);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write germline blacklist file: {}", e.toString());
        }
    }

    private void writeGermlineHotspots(final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes)
    {
        String germlineHotspotFile = formVersionFile(mOutputDir, "KnownHotspots.germline.vcf.gz", refGenomeVersion);

        String sageRefDir = mResourceRepoDir + SAGE_RESOURCE_DIR + File.separator + refGenomeVersion.identifier();
        String clinvarFile = formVersionFile(sageRefDir, "clinvar.vcf.gz", refGenomeVersion);

        GU_LOGGER.info("located clinvar file for {} at {}", refGenomeVersion, clinvarFile);

        List<String> germlineHotspotGenes = driverGenes.stream()
                .filter(x -> x.reportGermlineHotspot() != DriverGeneGermlineReporting.NONE)
                .map(x -> x.gene()).collect(Collectors.toList());

        GU_LOGGER.info("writing {} germline hotspots to {}", germlineHotspotGenes.size(), germlineHotspotFile);

        try
        {
            GermlineHotspotVCF.write(refGenomeVersion, clinvarFile, germlineHotspotFile, germlineHotspotGenes);

            GU_LOGGER.debug("germline hotspots VCF written");
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write germline hotspots file: {}", e.toString());
        }
    }

    public static void main(String[] args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addGenePanelOption(configBuilder, false);
        configBuilder.addPath(RESOURCE_REPO_DIR, true, RESOURCE_REPO_DIR_DESC);
        configBuilder.addConfigItem(PANEL_GENE_OVERRIDES, "List of comma-separated genes to include in panel");
        configBuilder.addConfigItem(GENE_ID_FILE, GENE_ID_FILE_DESC);
        addRefGenomeVersion(configBuilder);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenerateDriverGeneFiles generator = new GenerateDriverGeneFiles(configBuilder);
        generator.run();
    }
}
