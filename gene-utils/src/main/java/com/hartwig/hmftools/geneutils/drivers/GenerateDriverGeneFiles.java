package com.hartwig.hmftools.geneutils.drivers;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.genepanel.GeneNameMapping37to38;
import com.hartwig.hmftools.common.genome.genepanel.HmfExonPanelBed;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
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

    private static final String DRIVER_GENE_PANEL_TSV = "driver_gene_panel";
    private static final String RESOURCE_REPO_DIR = "resource_repo_dir";

    public GenerateDriverGeneFiles(final CommandLine cmd)
    {
        GU_LOGGER.info("starting driver gene panel generation");

        mDriverGenePanelFile = cmd.getOptionValue(DRIVER_GENE_PANEL_TSV);
        mResourceRepoDir = cmd.getOptionValue(RESOURCE_REPO_DIR);
        mOutputDir = parseOutputDir(cmd);
    }

    public void run() throws IOException
    {
        List<DriverGene> driverGenes = DriverGeneFile.read(mDriverGenePanelFile);
        GU_LOGGER.info("loaded {} driver genes from {}", driverGenes.size(), mDriverGenePanelFile);

        /*
        GeneNameMapping37to38 geneNameMapping = GeneNameMapping37to38.loadFromEmbeddedResource();
        List<DriverGene> v38DriverGenes = Lists.newArrayList();
        for(DriverGene input : driverGenes)
        {
            String v37Gene = input.gene();
            String v38Gene = geneNameMapping.v38Gene(input.gene());
            if(v37Gene.equals("LINC00290") || v37Gene.equals("LINC01001"))
            {
                v38DriverGenes.add(input);
            }
            else if(v38Gene.equals("NA"))
            {
                GU_LOGGER.debug(" Excluding: {}", v37Gene);
            }
            else
            {
                DriverGene converted = ImmutableDriverGene.builder().from(input).gene(v38Gene).build();
                v38DriverGenes.add(converted);
            }
        }
        */

        process(RefGenomeVersion.V37, driverGenes);
        process(RefGenomeVersion.V38, driverGenes);
    }

    public void process(final RefGenomeVersion refGenomeVersion, final List<DriverGene> driverGenes) throws IOException
    {
        String genePanelDir = mResourceRepoDir + "/gene_panel/" + (refGenomeVersion.is37() ? "37" : "38");
        String sageDir = mResourceRepoDir + "/sage/" + (refGenomeVersion.is37() ? "37" : "38");

        String qualityBedFile = getResourceURL(refGenomeVersion.addVersionToFilePath("/drivers/QualityRecalibration.bed"));
        Collection<GenomeRegion> qualityRecalibrationRegions = BEDFileLoader.fromBedFile(qualityBedFile).values();

        String clinvarFile = refGenomeVersion == RefGenomeVersion.V37
                ? mResourceRepoDir + "/sage/37/clinvar.37.vcf.gz"
                : mResourceRepoDir + "/sage/38/clinvar.38.vcf.gz";
        GU_LOGGER.info(" Located clinvar file for {} at {}", refGenomeVersion, clinvarFile);

        String driverGeneFile = refGenomeVersion.addVersionToFilePath(genePanelDir + "/DriverGenePanel.tsv");
        String somaticCodingWithoutUtr = refGenomeVersion.addVersionToFilePath(sageDir + "/ActionableCodingPanel.somatic.bed.gz");
        String germlineCodingWithUtr = refGenomeVersion.addVersionToFilePath(sageDir + "/ActionableCodingPanel.germline.bed.gz");
        String germlineCodingWithoutUtr = refGenomeVersion.addVersionToFilePath(sageDir + "/CoverageCodingPanel.germline.bed.gz");
        String germlineHotspotFile = refGenomeVersion.addVersionToFilePath(sageDir + "/KnownHotspots.germline.vcf.gz");
        String germlineSliceFile = refGenomeVersion.addVersionToFilePath(sageDir + "/SlicePanel.germline.bed.gz");
        String germlineBlacklistFile = refGenomeVersion.addVersionToFilePath(sageDir + "/KnownBlacklist.germline.vcf.gz");

        Collections.sort(driverGenes);

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
        GU_LOGGER.info("writing {} {} bed file at {}", refGenomeVersion, "somatic", somaticCodingWithoutUtr);
        createNamedBedFiles(false, somaticCodingWithoutUtr, somaticGenes, transcripts);

        GU_LOGGER.info("writing {} {} coverage bed file at {}", refGenomeVersion, "germline", germlineCodingWithoutUtr);
        createNamedBedFiles(false, germlineCodingWithoutUtr, allGenes, transcripts);

        GU_LOGGER.info("writing {} {} bed file at {}", refGenomeVersion, "germline", germlineCodingWithUtr);
        List<GenomeRegion> germlinePanel = createUnnamedBedFiles(true, germlineCodingWithUtr, allGenes, transcripts);

        GU_LOGGER.info("writing {} germline hotspots at {}", refGenomeVersion, germlineHotspotFile);
        List<GenomeRegion> germlineHotspots =
                new GermlineHotspotVCF(germlineHotspotGenes(driverGenes)).process(clinvarFile, germlineHotspotFile);

        GU_LOGGER.info("writing {} germline slice file at {}", refGenomeVersion, germlineSliceFile);
        createSliceFile(germlineSliceFile, germlinePanel, germlineHotspots, qualityRecalibrationRegions);

        GU_LOGGER.info("writing {} germline blacklist file at {}", refGenomeVersion, germlineBlacklistFile);
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

    private static void createNamedBedFiles(
            boolean includeUTR, final String file, final Set<String> genes,
            final List<HmfTranscriptRegion> transcripts) throws IOException
    {
        List<NamedBed> somaticBed = HmfExonPanelBed.createNamedCodingRegions(includeUTR, genes, transcripts);
        NamedBedFile.writeBedFile(file, somaticBed);
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
