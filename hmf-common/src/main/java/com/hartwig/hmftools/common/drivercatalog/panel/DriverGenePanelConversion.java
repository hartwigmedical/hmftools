package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedFile;
import com.hartwig.hmftools.common.genome.genepanel.HmfExonPanelBed;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class DriverGenePanelConversion {

    private static final String OUTPUT_DIR = "/Users/jon/hmf/resources/newGenePanel/";
    private static final String INPUT_DIR = "/Users/jon/hmf/resources/";
    private static final Logger LOGGER = LogManager.getLogger(DriverGenePanelConversion.class);

    public static void main(String[] args) throws IOException {
        LOGGER.info("Starting driver gene panel generation");

        String templateFile = OUTPUT_DIR + "DriverGenePanel.template.txt";

        DndsGeneNameMap geneNameMap = new DndsGeneNameMap();
        List<DriverGene> inputDriverGenes = DriverGeneFile.read(templateFile);
        List<DriverGene> outputDriverGenes = Lists.newArrayList();
        for (DriverGene input : inputDriverGenes) {
            final String hg19Gene = input.gene();
            final String hg38Gene = geneNameMap.hg38Gene(input.gene());
            if (hg19Gene.equals("LINC00290") || hg19Gene.equals("LINC01001")) {
                outputDriverGenes.add(input);
            } else if (hg38Gene.equals("NA")) {
                System.out.println("Excluding: " + hg19Gene);
            } else {
                DriverGene converted = ImmutableDriverGene.builder().from(input).gene(hg38Gene).build();
                outputDriverGenes.add(converted);
            }
        }

        process(DriverGenePanelAssembly.HG19, inputDriverGenes);
        process(DriverGenePanelAssembly.HG38, outputDriverGenes);

    }

    private static void process(@NotNull final DriverGenePanelAssembly assembly, @NotNull final List<DriverGene> driverGenes)
            throws IOException {
        final String qualityBedFile = assembly == DriverGenePanelAssembly.HG19
                ? getResourceURL("/drivercatalog/QualityRecalibration.hg19.bed")
                : getResourceURL("/drivercatalog/QualityRecalibration.hg38.bed");
        Collection<GenomeRegion> qualityRecalibrationRegions = BEDFileLoader.fromBedFile(qualityBedFile).values();

        final String extension = assembly.toString().toLowerCase();
        final String driverGeneFile = String.format("%s/DriverGenePanel.%s.tsv", OUTPUT_DIR, extension);
        final String somaticCodingWithoutUtr = String.format("%s/ActionableCodingPanel.somatic.%s.bed", OUTPUT_DIR, extension);
        final String germlineCodingWithUtr = String.format("%s/ActionableCodingPanel.germline.%s.bed", OUTPUT_DIR, extension);
        final String germlineCodingWithoutUtr = String.format("%s/CoverageCodingPanel.germline.%s.bed", OUTPUT_DIR, extension);
        final String germlineHotspotFile = String.format("%s/KnownHotspots.germline.%s.vcf.gz", OUTPUT_DIR, extension);
        final String germlineSliceFile = String.format("%s/SlicePanel.germline.%s.bed", OUTPUT_DIR, extension);
        final String clinvarFile = String.format("%s/clinvar.%s.vcf.gz", INPUT_DIR, extension);
        final String germlineBlacklistFile = String.format("%s/KnownBlacklist.germline.%s.vcf.gz", OUTPUT_DIR, extension);

        Collections.sort(driverGenes);

        // Validate
        DriverGenePanelFactory.create(assembly, driverGenes);

        // Write out driver gene panel
        DriverGeneFile.write(driverGeneFile, driverGenes);

        final Set<String> germlineGenes = germlineGenes(driverGenes);
        final Set<String> somaticGenes = somaticGenes(driverGenes);
        final Set<String> allGenes = Sets.newHashSet();
        allGenes.addAll(germlineGenes);
        allGenes.addAll(somaticGenes);

        // Write out actionable bed files
        final List<HmfTranscriptRegion> transcripts = assembly.equals(DriverGenePanelAssembly.HG19)
                ? HmfGenePanelSupplier.allGeneList37()
                : HmfGenePanelSupplier.allGeneList38();

        // Write out driver gene panel
        LOGGER.info("Creating {} {} bed file", assembly, "somatic");
        createNamedBedFiles(false, somaticCodingWithoutUtr, somaticGenes, transcripts);

        LOGGER.info("Creating {} {} coverage bed file", assembly, "germline");
        createNamedBedFiles(false, germlineCodingWithoutUtr, allGenes, transcripts);

        LOGGER.info("Creating {} {} bed file", assembly, "germline");
        List<GenomeRegion> germlinePanel = createUnnamedBedFiles(true, germlineCodingWithUtr, allGenes, transcripts);

        // Write out germline hotspot files
        List<GenomeRegion> germlineHotspots =
                new GermlineHotspotVCF(germlineHotspotGenes(driverGenes)).process(clinvarFile, germlineHotspotFile);

        // Write out germline slice file.
        createSliceFile(germlineSliceFile, germlinePanel, germlineHotspots, qualityRecalibrationRegions);

        // Write out germline blacklist file.
        final List<VariantContext> germlineBlackList = assembly == DriverGenePanelAssembly.HG19
                ? GermlineResources.grch37Blacklist()
                : GermlineResources.grch38Blacklist();
        GermlineBlacklistVCF.process(germlineBlacklistFile, germlineBlackList );

    }

    private static void createSliceFile(String file, Collection<GenomeRegion>... allRegions) throws IOException {
        GenomeRegionsBuilder builder = new GenomeRegionsBuilder();
        for (Collection<GenomeRegion> regions : allRegions) {
            for (GenomeRegion region : regions) {
                builder.addRegion(region.chromosome(), region.start() - 100, region.end() + 100);
            }
        }

        List<GenomeRegion> combined = builder.build();
        NamedBedFile.writeUnnamedBedFile(file, combined);
    }

    private static void createNamedBedFiles(boolean includeUTR, String file, Set<String> genes, List<HmfTranscriptRegion> transcripts)
            throws IOException {
        final List<NamedBed> somaticBed = HmfExonPanelBed.createNamedCodingRegions(includeUTR, genes, transcripts);
        NamedBedFile.writeBedFile(file, somaticBed);
    }

    private static List<GenomeRegion> createUnnamedBedFiles(boolean includeUTR, String file, Set<String> genes,
            List<HmfTranscriptRegion> transcripts) throws IOException {
        final List<GenomeRegion> somaticBed = HmfExonPanelBed.createUnnamedCodingRegions(includeUTR, genes, transcripts);
        NamedBedFile.writeUnnamedBedFile(file, somaticBed);
        return somaticBed;
    }

    @NotNull
    static Set<String> somaticGenes(@NotNull final List<DriverGene> genePanel) {
        final Set<String> actionableGenes = Sets.newHashSet();
        for (DriverGene driverGene : genePanel) {
            if (driverGene.reportSomatic()) {
                actionableGenes.add(driverGene.gene());
            }
        }

        return actionableGenes;
    }

    @NotNull
    static Set<String> germlineGenes(@NotNull final List<DriverGene> genePanel) {
        final Set<String> actionableGenes = Sets.newHashSet();
        for (DriverGene driverGene : genePanel) {
            if (driverGene.reportGermline()) {
                actionableGenes.add(driverGene.gene());
            }

        }

        return actionableGenes;
    }

    @NotNull
    static Set<String> germlineHotspotGenes(@NotNull final List<DriverGene> genePanel) {
        final Set<String> actionableGenes = Sets.newHashSet();
        for (DriverGene driverGene : genePanel) {
            if (driverGene.reportGermlineHotspot() != DriverGeneGermlineReporting.NONE) {
                actionableGenes.add(driverGene.gene());
            }

        }

        return actionableGenes;
    }

    private static String getResourceURL(String location) {
        return DriverGenePanelConversion.class.getResource(location).toString();
    }

}
