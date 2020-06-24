package com.hartwig.hmftools.serve.vicc;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.util.ProteinKeyFormatter;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.curation.CurationKey;
import com.hartwig.hmftools.serve.vicc.curation.FeatureCurator;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeExtractor;
import com.hartwig.hmftools.serve.vicc.signatures.SignaturesExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
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

    private static final Integer MAX_ENTRIES = null;

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String viccJsonPath;
        String refGenomeFastaFile;
        boolean generateHotspots;
        String hotspotVcf;

        if (hostname.toLowerCase().contains("datastore")) {
            viccJsonPath = "/data/common/dbs/vicc/all.json";
            refGenomeFastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            generateHotspots = true;
            hotspotVcf = System.getProperty("user.home") + "/tmp/hotspotsVicc.vcf";
        } else {
            viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";
            refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            generateHotspots = false;
            hotspotVcf = System.getProperty("user.home") + "/hmf/tmp/hotspotsVicc.vcf";
        }

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        LOGGER.debug("Configured '{}' as the VICC json path", viccJsonPath);
        LOGGER.debug("Configured '{}' as the reference fasta path", refGenomeFastaFile);
        LOGGER.debug("Configured '{}' as the hotspot output VCF", hotspotVcf);
        LOGGER.debug("Configured '{}' for generating hotspots yes/no", generateHotspots);

        List<ViccSource> sources = Lists.newArrayList(ViccSource.ONCOKB);
        ViccQuerySelection querySelection =
                ImmutableViccQuerySelection.builder().sourcesToFilterOn(sources).maxEntriesToInclude(MAX_ENTRIES).build();
        List<ViccEntry> viccEntries = readAndCurate(viccJsonPath, querySelection);

        HotspotExtractor hotspotExtractor =
                generateHotspots ? HotspotExtractor.transvarWithRefGenome(refGenomeVersion, refGenomeFastaFile) : HotspotExtractor.dummy();

        ViccExtractor viccExtractor = new ViccExtractor(hotspotExtractor,
                new CopyNumberExtractor(),
                new FusionExtractor(),
                new GeneLevelEventExtractor(),
                new GeneRangeExtractor(),
                new SignaturesExtractor());

        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = viccExtractor.extractFromViccEntries(viccEntries);

        analyzeExtractionResults(resultsPerEntry);

        if (generateHotspots) {
            writeHotspots(hotspotVcf, resultsPerEntry);
        }
    }

    @NotNull
    private static List<ViccEntry> readAndCurate(@NotNull String viccJsonPath, @NotNull ViccQuerySelection querySelection)
            throws IOException {
        FeatureCurator curator = new FeatureCurator();

        LOGGER.info("Reading VICC json from '{}' with sources '{}'", viccJsonPath, querySelection.sourcesToFilterOn());
        List<ViccEntry> viccEntries = curate(ViccJsonReader.readSelection(viccJsonPath, querySelection), curator);
        LOGGER.info(" Read and curated {} entries", viccEntries.size());

        LOGGER.info("Analyzing usage of curation configuration keys");
        int issueCount = 0;
        for (Map.Entry<ViccSource, Set<CurationKey>> entry : curator.unusedCurationKeysPerSource().entrySet()) {
            ViccSource source = entry.getKey();
            Set<CurationKey> unusedKeys = entry.getValue();
            if (!unusedKeys.isEmpty()) {
                LOGGER.warn("Found {} unused curation configuration entries for {}", unusedKeys.size(), source);
                for (CurationKey unusedKey : unusedKeys) {
                    issueCount++;
                    LOGGER.warn(" - {}", unusedKey);
                }
            }
        }
        LOGGER.info("Finished analyzing usage of curation configuration keys. Found {} issues", issueCount);

        return viccEntries;
    }

    @NotNull
    private static List<ViccEntry> curate(@NotNull List<ViccEntry> viccEntries, @NotNull FeatureCurator curator) {
        // TODO: Move to a function inside the curator and remove entries which have a blacklisted feature.
        List<ViccEntry> curatedViccEntries = Lists.newArrayList();

        for (ViccEntry entry : viccEntries) {
            ImmutableViccEntry.Builder builder = ImmutableViccEntry.builder().from(entry);
            List<Feature> curatedFeatures = Lists.newArrayList();
            for (Feature feature : entry.features()) {
                Feature curatedFeature = curator.curate(entry, feature);
                if (curatedFeature != null) {
                    curatedFeatures.add(curator.curate(entry, feature));
                }
            }
            curatedViccEntries.add(builder.features(curatedFeatures).build());
        }
        return curatedViccEntries;
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

        LOGGER.info("No genomic events derived for {} features", featuresWithoutGenomicEvents.size());
        for (Feature feature : featuresWithoutGenomicEvents) {
            LOGGER.debug(" No genomic events derived from '{}' in '{}'", feature.name(), feature.geneSymbol());
        }
    }

    private static void writeHotspots(@NotNull String hotspotVcf, @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(hotspotVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        writer.writeHeader(header);

        for (Map.Entry<VariantHotspot, String> entry : convertAndSort(resultsPerEntry).entrySet()) {
            VariantHotspot hotspot = entry.getKey();
            List<Allele> hotspotAlleles = buildAlleles(hotspot);

            VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                    .source("VICC")
                    .chr(hotspot.chromosome())
                    .start(hotspot.position())
                    .alleles(hotspotAlleles)
                    .computeEndFromAlleles(hotspotAlleles, (int) hotspot.position())
                    .attribute("feature", entry.getValue())
                    .make();

            LOGGER.debug("Writing {}", variantContext);
            writer.add(variantContext);

        }
        writer.close();
    }

    @NotNull
    private static List<Allele> buildAlleles(@NotNull VariantHotspot hotspot) {
        Allele ref = Allele.create(hotspot.ref(), true);
        Allele alt = Allele.create(hotspot.alt(), false);

        return Lists.newArrayList(ref, alt);
    }

    @NotNull
    private static Map<VariantHotspot, String> convertAndSort(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        Map<VariantHotspot, String> convertedMap = Maps.newTreeMap(new VariantHotspotComparator());
        for (Map.Entry<ViccEntry, ViccExtractionResult> entryResult : resultsPerEntry.entrySet()) {
            for (Map.Entry<Feature, List<VariantHotspot>> featureResult : entryResult.getValue().hotspotsPerFeature().entrySet()) {
                String newFeatureAttribute = toAttribute(entryResult.getKey(), featureResult.getKey());
                for (VariantHotspot hotspot : featureResult.getValue()) {
                    String existingFeatureAttribute = convertedMap.get(hotspot);
                    if (existingFeatureAttribute != null && !existingFeatureAttribute.equals(newFeatureAttribute)) {
                        if (!existingFeatureAttribute.contains("null") && !newFeatureAttribute.contains("null")) {
                            LOGGER.warn("Hotspot {} already recorded but for different feature! Existing feature={}, new feature={}",
                                    hotspot,
                                    existingFeatureAttribute,
                                    newFeatureAttribute);
                        }
                    } else if (existingFeatureAttribute == null || existingFeatureAttribute.contains("null")) {
                        // We favor feature attributes which are based on explicit transcript.
                        convertedMap.put(hotspot, newFeatureAttribute);
                    }
                }
            }
        }

        return convertedMap;
    }

    @NotNull
    private static String toAttribute(@NotNull ViccEntry viccEntry, @NotNull Feature feature) {
        return ProteinKeyFormatter.toProteinKey(feature.geneSymbol(), viccEntry.transcriptId(), feature.proteinAnnotation());
    }
}
