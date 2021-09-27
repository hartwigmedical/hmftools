package com.hartwig.hmftools.geneutils.drivers;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.genepanel.HmfExonPanelBed;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.genepanel.GeneNameMapping;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.variant.variantcontext.VariantContext;

public class GenerateDriverGeneFiles
{
    private static final Logger LOGGER = LogManager.getLogger(GenerateDriverGeneFiles.class);

    private static final String NEW_DRIVER_GENE_PANEL_37_TSV = "new_driver_gene_panel_37_tsv";
    private static final String RESOURCE_REPO_DIR = "resource_repo_dir";

    public static void main(String[] args) throws IOException, ParseException
    {
        // See also https://github.com/hartwigmedical/scratchpad/tree/master/genePanel
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        LOGGER.info("Starting driver gene panel generation");

        new GenerateDriverGeneFiles(cmd.getOptionValue(NEW_DRIVER_GENE_PANEL_37_TSV), cmd.getOptionValue(RESOURCE_REPO_DIR)).run();

        LOGGER.info("Complete");
    }

    private final String mNewDriverGenePanel37Tsv;
    private final String mResourceRepoDir;

    public GenerateDriverGeneFiles(final String newDriverGenePanel37Tsv, final String resourceRepoDir)
    {
        mNewDriverGenePanel37Tsv = newDriverGenePanel37Tsv;
        mResourceRepoDir = resourceRepoDir;
    }

    public void run() throws IOException
    {
        List<DriverGene> v37DriverGenes = DriverGeneFile.read(mNewDriverGenePanel37Tsv);
        LOGGER.info(" Loaded {} driver genes from {}", v37DriverGenes.size(), mNewDriverGenePanel37Tsv);

        GeneNameMapping geneNameMapping = GeneNameMapping.loadFromEmbeddedResource();
        List<DriverGene> v38DriverGenes = Lists.newArrayList();
        for(DriverGene input : v37DriverGenes)
        {
            String v37Gene = input.gene();
            String v38Gene = geneNameMapping.v38Gene(input.gene());
            if(v37Gene.equals("LINC00290") || v37Gene.equals("LINC01001"))
            {
                v38DriverGenes.add(input);
            }
            else if(v38Gene.equals("NA"))
            {
                LOGGER.debug(" Excluding: {}", v37Gene);
            }
            else
            {
                DriverGene converted = ImmutableDriverGene.builder().from(input).gene(v38Gene).build();
                v38DriverGenes.add(converted);
            }
        }

        process(RefGenomeVersion.V37, v37DriverGenes);
        process(RefGenomeVersion.V38, v38DriverGenes);
    }

    public void process(final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes) throws IOException
    {
        String genePanelDir = mResourceRepoDir + "/gene_panel/" + (refGenomeVersion.is37() ? "37" : "38");
        String sageDir = mResourceRepoDir + "/sage/" + (refGenomeVersion.is37() ? "37" : "38");

        String qualityBedFile = getResourceURL(refGenomeVersion.addVersionToFilePath("/drivercatalog/QualityRecalibration.bed"));
        Collection<GenomeRegion> qualityRecalibrationRegions = BEDFileLoader.fromBedFile(qualityBedFile).values();

        String clinvarFile = refGenomeVersion == RefGenomeVersion.V37
                ? mResourceRepoDir + "/sage/37/clinvar.37.vcf.gz"
                : mResourceRepoDir + "/sage/38/clinvar.38.vcf.gz";
        LOGGER.info(" Located clinvar file for {} at {}", refGenomeVersion, clinvarFile);

        String driverGeneFile = refGenomeVersion.addVersionToFilePath(genePanelDir + "/DriverGenePanel.tsv");
        String somaticCodingWithoutUtr = refGenomeVersion.addVersionToFilePath(sageDir + "/ActionableCodingPanel.somatic.bed.gz");
        String germlineCodingWithUtr = refGenomeVersion.addVersionToFilePath(sageDir + "/ActionableCodingPanel.germline.bed.gz");
        String germlineCodingWithoutUtr = refGenomeVersion.addVersionToFilePath(sageDir + "/CoverageCodingPanel.germline.bed.gz");
        String germlineHotspotFile = refGenomeVersion.addVersionToFilePath(sageDir + "/KnownHotspots.germline.vcf.gz");
        String germlineSliceFile = refGenomeVersion.addVersionToFilePath(sageDir + "/SlicePanel.germline.bed.gz");
        String germlineBlacklistFile = refGenomeVersion.addVersionToFilePath(sageDir + "/KnownBlacklist.germline.vcf.gz");

        Collections.sort(driverGenes);

        // This will throw an exception if there is a problem with the driver genes for this ref genome version.
        // DriverGenePanelFactory.create(refGenomeVersion, driverGenes);

        // Write out driver gene panel
        DriverGeneFile.write(driverGeneFile, driverGenes);

        Set<String> germlineGenes = germlineGenes(driverGenes);
        Set<String> somaticGenes = somaticGenes(driverGenes);
        Set<String> allGenes = Sets.newHashSet();
        allGenes.addAll(germlineGenes);
        allGenes.addAll(somaticGenes);

        // Write out actionable bed files
        List<HmfTranscriptRegion> transcripts =
                refGenomeVersion == RefGenomeVersion.V37 ? HmfGenePanelSupplier.allGeneList37() : HmfGenePanelSupplier.allGeneList38();

        // Write out driver gene panel
        LOGGER.info(" Writing {} {} bed file at {}", refGenomeVersion, "somatic", somaticCodingWithoutUtr);
        createNamedBedFiles(false, somaticCodingWithoutUtr, somaticGenes, transcripts);

        LOGGER.info(" Writing {} {} coverage bed file at {}", refGenomeVersion, "germline", germlineCodingWithoutUtr);
        createNamedBedFiles(false, germlineCodingWithoutUtr, allGenes, transcripts);

        LOGGER.info(" Writing {} {} bed file at {}", refGenomeVersion, "germline", germlineCodingWithUtr);
        List<GenomeRegion> germlinePanel = createUnnamedBedFiles(true, germlineCodingWithUtr, allGenes, transcripts);

        LOGGER.info(" Writing {} germline hotspots at {}", refGenomeVersion, germlineHotspotFile);
        List<GenomeRegion> germlineHotspots =
                new GermlineHotspotVCF(germlineHotspotGenes(driverGenes)).process(clinvarFile, germlineHotspotFile);

        LOGGER.info(" Writing {} germline slice file at {}", refGenomeVersion, germlineSliceFile);
        createSliceFile(germlineSliceFile, germlinePanel, germlineHotspots, qualityRecalibrationRegions);

        LOGGER.info(" Writing {} germline blacklist file at {}", refGenomeVersion, germlineBlacklistFile);
        List<VariantContext> germlineBlackList =
                refGenomeVersion == RefGenomeVersion.V37 ? GermlineResources.blacklist37() : GermlineResources.blacklist38();
        GermlineBlacklistVCF.process(germlineBlacklistFile, germlineBlackList);
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

    private static void createNamedBedFiles(boolean includeUTR, final String file, final Set<String> genes,
            final List<HmfTranscriptRegion> transcripts) throws IOException
    {
        List<NamedBed> somaticBed = HmfExonPanelBed.createNamedCodingRegions(includeUTR, genes, transcripts);
        NamedBedFile.writeBedFile(file, somaticBed);
    }

    private static List<GenomeRegion> createUnnamedBedFiles(boolean includeUTR, final String file, final Set<String> genes,
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

    private static Options createOptions()
    {
        Options options = new Options();

        options.addOption(NEW_DRIVER_GENE_PANEL_37_TSV, true, "File containing the driver gene panel for 37");
        options.addOption(RESOURCE_REPO_DIR, true, "The directory holding the public hmf resources repo");

        return options;
    }
}
