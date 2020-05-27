package com.hartwig.hmftools.serve.vicc;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeExtractor;
import com.hartwig.hmftools.serve.vicc.signatures.SignaturesExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;
import com.hartwig.hmftools.vicc.selection.ImmutableViccQuerySelection;
import com.hartwig.hmftools.vicc.selection.ViccQuerySelection;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class ViccExtractorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccExtractorTestApplication.class);

    private static final boolean RUN_ON_SERVER = false;
    private static final boolean TRANSVAR_ENABLED = true;
    private static final boolean WRITE_HOTSPOTS_TO_VCF = true;

    public static void main(String[] args) throws IOException, InterruptedException {
        Configurator.setRootLevel(Level.DEBUG);

        String viccJsonPath;
        String refGenomeFastaFile;
        String hotspotVcf;

        if (RUN_ON_SERVER) {
            viccJsonPath = "/data/common/dbs/vicc/all.json";
            refGenomeFastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            hotspotVcf = System.getProperty("user.home") + "/hotspotsVicc.vcf";
        } else {
            viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";
            refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            hotspotVcf = System.getProperty("user.home") + "/hmf/tmp/hotspotsVicc.vcf";
        }

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;

        ViccSource source = ViccSource.ONCOKB;
        LOGGER.info("Reading VICC json from {} with source '{}'", viccJsonPath, source);
        ViccQuerySelection querySelection = ImmutableViccQuerySelection.builder().addSourcesToFilterOn(source).maxEntriesToInclude(2).
                build();
        List<ViccEntry> viccEntries = ViccJsonReader.readSelection(viccJsonPath, querySelection);
        LOGGER.info(" Read {} entries", viccEntries.size());

        ViccExtractor viccExtractor =
                new ViccExtractor(HotspotExtractor.withRefGenome(refGenomeVersion, refGenomeFastaFile, TRANSVAR_ENABLED),
                        new CopyNumberExtractor(),
                        new FusionExtractor(),
                        new GeneLevelEventExtractor(),
                        new GeneRangeExtractor(),
                        new SignaturesExtractor());

        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = viccExtractor.extractFromViccEntries(viccEntries);

        analyzeExtractionResults(resultsPerEntry);

        if (WRITE_HOTSPOTS_TO_VCF) {
            writeHotspots(hotspotVcf, resultsPerEntry.values());
        }
    }

    private static void analyzeExtractionResults(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        List<Feature> featuresWithoutGenomicEvents = Lists.newArrayList();
        int totalFeatureCount = 0;
        int featuresWithHotspotsCount = 0;
        int totalHotspotsCount = 0;
        int featuresWithCopyNumberCount = 0;
        int featuresWithFusionCount = 0;
        int featuresWithGeneLevelEventCount = 0;
        int featuresWithGeneRangeCount = 0;
        int featuresWithSignatureCount = 0;

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            ViccEntry viccEntry = entry.getKey();
            ViccExtractionResult viccExtractionResult = entry.getValue();
            for (Feature feature : viccEntry.features()) {
                List<VariantHotspot> hotspotsForFeature = viccExtractionResult.hotspotsPerFeature().get(feature);
                KnownAmplificationDeletion ampDelForFeature = viccExtractionResult.ampsDelsPerFeature().get(feature);
                String fusionForFeature = viccExtractionResult.fusionsPerFeature().get(feature);
                String geneLevelEventForFeature = viccExtractionResult.geneLevelEventsPerFeature().get(feature);
                String geneRangeForFeature = viccExtractionResult.geneRangesPerFeature().get(feature);
                String signatureForFeature = viccExtractionResult.signaturesPerFeature().get(feature);

                if (hotspotsForFeature == null && ampDelForFeature == null && fusionForFeature == null && geneLevelEventForFeature == null
                        && geneRangeForFeature == null && signatureForFeature == null) {
                    featuresWithoutGenomicEvents.add(feature);
                } else {
                    if (hotspotsForFeature != null) {
                        featuresWithHotspotsCount++;
                        totalHotspotsCount += hotspotsForFeature.size();
                    }

                    if (ampDelForFeature != null) {
                        featuresWithCopyNumberCount++;
                    }

                    if (fusionForFeature != null) {
                        featuresWithFusionCount++;
                    }

                    if (geneLevelEventForFeature != null) {
                        featuresWithGeneLevelEventCount++;
                    }

                    if (geneRangeForFeature != null) {
                        featuresWithGeneRangeCount++;
                    }

                    if (signatureForFeature != null) {
                        featuresWithSignatureCount++;
                    }
                }

                totalFeatureCount++;
            }
        }
        LOGGER.info("Extraction performed on {} features from {} entries", totalFeatureCount, resultsPerEntry.size());
        LOGGER.info(" Extracted {} hotspots for {} features", totalHotspotsCount, featuresWithHotspotsCount);
        LOGGER.info(" Extracted {} known amps and dels", featuresWithCopyNumberCount);
        LOGGER.info(" Extracted {} fusions", featuresWithFusionCount);
        LOGGER.info(" Extracted {} gene level events", featuresWithGeneLevelEventCount);
        LOGGER.info(" Extracted {} gene ranges", featuresWithGeneRangeCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);

        LOGGER.info("No genomic events found for {} features", featuresWithoutGenomicEvents.size());
        for (Feature feature : featuresWithoutGenomicEvents) {
            LOGGER.debug(" {} in {}: {}", feature.name(), feature.geneSymbol(), feature);
        }
    }

    private static void writeHotspots(@NotNull String hotspotVcf, @NotNull Collection<ViccExtractionResult> extractionResults) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(hotspotVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        writer.writeHeader(header);

        for (ViccExtractionResult result : extractionResults) {
            for (Map.Entry<Feature, List<VariantHotspot>> entry : result.hotspotsPerFeature().entrySet()) {
                String featureAttribute = entry.getKey().geneSymbol() + ":p." + entry.getKey().name();

                for (VariantHotspot hotspot : entry.getValue()) {
                    List<Allele> hotspotAlleles = buildAlleles(hotspot);

                    VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                            .source("VICC")
                            .chr(hotspot.chromosome())
                            .start(hotspot.position())
                            .alleles(hotspotAlleles)
                            .computeEndFromAlleles(hotspotAlleles, (int) hotspot.position())
                            .attribute("feature", featureAttribute)
                            .make();

                    LOGGER.debug("Writing {}", variantContext);
                    writer.add(variantContext);
                }
            }
        }
        writer.close();
    }

    @NotNull
    private static List<Allele> buildAlleles(@NotNull VariantHotspot hotspot) {
        Allele ref = Allele.create(hotspot.ref(), true);
        Allele alt = Allele.create(hotspot.alt(), false);

        return Lists.newArrayList(ref, alt);
    }

}
