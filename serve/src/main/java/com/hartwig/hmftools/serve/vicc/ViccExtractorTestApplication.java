package com.hartwig.hmftools.serve.vicc;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.net.InetAddress;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.util.ProteinKeyFormatter;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.curation.CurationKey;
import com.hartwig.hmftools.serve.vicc.curation.FeatureCurator;
import com.hartwig.hmftools.serve.vicc.fusion.FusionAnnotation;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeAnnotation;
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

    private static final Integer MAX_ENTRIES = null;

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String viccJsonPath;
        String refGenomeFastaFile;
        boolean generateHotspots;
        String hotspotVcf = null;
        String rangesTsv = null;
        String fusionTsv = null;


        if (hostname.toLowerCase().contains("datastore")) {
            viccJsonPath = "/data/common/dbs/vicc/all.json";
            refGenomeFastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            generateHotspots = false; //for generate hostpots set to true
            hotspotVcf = System.getProperty("user.home") + "/tmp/hotspotsVicc.vcf";
            rangesTsv = System.getProperty("user.home") + "/tmp/rangesVicc.vcf";
            fusionTsv = System.getProperty("user.home") + "/tmp/fusionVicc.vcf";


        } else {
            viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";
            refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            generateHotspots = false;
        }

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        LOGGER.debug("Configured '{}' as the VICC json path", viccJsonPath);
        LOGGER.debug("Configured '{}' as the reference fasta path", refGenomeFastaFile);
        LOGGER.debug("Configured '{}' as the hotspot output VCF", hotspotVcf);
        LOGGER.debug("Configured '{}' as the ranges output TSV", rangesTsv);
        LOGGER.debug("Configured '{}' as the fusion output TSV", fusionTsv);
        LOGGER.debug("Configured '{}' for generating hotspots yes/no", generateHotspots);

        List<ViccSource> sources = Lists.newArrayList(ViccSource.CIVIC, ViccSource.JAX, ViccSource.ONCOKB, ViccSource.CGI);
        ViccQuerySelection querySelection =
                ImmutableViccQuerySelection.builder().sourcesToFilterOn(sources).maxEntriesToInclude(MAX_ENTRIES).build();
        List<ViccEntry> viccEntries = readAndCurate(viccJsonPath, querySelection);

        HotspotExtractor hotspotExtractor =
                generateHotspots ? HotspotExtractor.transvarWithRefGenome(refGenomeVersion, refGenomeFastaFile) : HotspotExtractor.dummy();

        Map<String, HmfTranscriptRegion> transcriptPerGeneMap = HmfGenePanelSupplier.allGenesMap37();

        ViccExtractor viccExtractor = new ViccExtractor(hotspotExtractor,
                new CopyNumberExtractor(),
                new FusionExtractor(),
                new GeneLevelEventExtractor(),
                new GeneRangeExtractor(transcriptPerGeneMap),
                new SignaturesExtractor());

        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = viccExtractor.extractFromViccEntries(viccEntries);

        analyzeExtractionResults(resultsPerEntry);

        if (generateHotspots && hotspotVcf != null) {
            writeHotspots(hotspotVcf, resultsPerEntry);
        }
        if (rangesTsv != null) {
            writeRanges(rangesTsv, resultsPerEntry);
        }
        if (fusionTsv != null) {
            writeFusion(fusionTsv, resultsPerEntry);
        }
    }

    @NotNull
    private static List<ViccEntry> readAndCurate(@NotNull String viccJsonPath, @NotNull ViccQuerySelection querySelection)
            throws IOException {
        FeatureCurator curator = new FeatureCurator();

        LOGGER.info("Reading VICC json from '{}' with sources '{}'", viccJsonPath, querySelection.sourcesToFilterOn());
        List<ViccEntry> viccEntries = ViccJsonReader.readSelection(viccJsonPath, querySelection);
        LOGGER.info(" Read {} entries. Starting curation", viccEntries.size());

        List<ViccEntry> curatedViccEntries = curator.curate(viccEntries);
        LOGGER.info(" Finished curation. {} curated entries remaining. {} entries have been removed.",
                curatedViccEntries.size(),
                viccEntries.size() - curatedViccEntries.size());

        LOGGER.info("Analyzing usage of curation configuration keys");
        Set<CurationKey> unusedCurationKeys = curator.unusedCurationKeys();
        for (CurationKey unusedKey : unusedCurationKeys) {
            LOGGER.warn(" Unused key found: '{}'", unusedKey);
        }

        LOGGER.info("Finished analyzing usage of curation configuration keys. Found {} unused configuration entries.",
                unusedCurationKeys.size());

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
                FusionAnnotation fusionForFeature = viccExtractionResult.fusionsPerFeature().get(feature);
                String geneLevelEventForFeature = viccExtractionResult.geneLevelEventsPerFeature().get(feature);
                GeneRangeAnnotation geneRangeForFeature = viccExtractionResult.geneRangesPerFeature().get(feature);
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
                        // LOGGER.debug("Feature '{}' in '{}' interpreted as amp/del", feature.name(), feature.geneSymbol());
                        featuresWithCopyNumberCount++;
                    }

                    if (fusionForFeature != null) {
//                        LOGGER.debug("Feature '{}' in '{}' interpreted as ''{}",
//                                fusionForFeature.fusion(),
//                                feature.geneSymbol(),
//                                fusionForFeature.fusionEvent());
                        featuresWithFusionCount++;
                    }

                    if (geneLevelEventForFeature != null) {
                        //  LOGGER.debug("Feature '{}' in '{}' interpreted as gene level event", feature.name(), feature.geneSymbol());
                        featuresWithGeneLevelEventCount++;
                    }

                    if (geneRangeForFeature != null) {
                        //  LOGGER.debug("Feature '{}' in '{}' interpreted as gene range event", feature.name(), feature.geneSymbol());
                        featuresWithGeneRangeCount++;
                    }

                    if (signatureForFeature != null) {
                         // LOGGER.debug("Feature '{}' in '{}' interpreted as signature event", feature.name(), feature.geneSymbol());
                        featuresWithSignatureCount++;
                    }
                }

                totalFeatureCount++;
            }
        }

       // LOGGER.info("No genomic events derived for {} features.", featuresWithoutGenomicEvents.size());
        for (Feature feature : featuresWithoutGenomicEvents) {
            if (!FeatureIgnoreUtil.canIgnore(feature)) {
//                LOGGER.debug(" No genomic events derived from '{}' in '{}'", feature.name(), feature.geneSymbol());
//                  LOGGER.info(feature);
            }
        }

        LOGGER.info("Extraction performed on {} features from {} entries", totalFeatureCount, resultsPerEntry.size());
        LOGGER.info(" Extracted {} hotspots for {} features", totalHotspotsCount, featuresWithHotspotsCount);
        LOGGER.info(" Extracted {} known amps and dels", featuresWithCopyNumberCount);
        LOGGER.info(" Extracted {} fusions", featuresWithFusionCount);
        LOGGER.info(" Extracted {} gene level events", featuresWithGeneLevelEventCount);
        LOGGER.info(" Extracted {} gene ranges", featuresWithGeneRangeCount);
        LOGGER.info(" Extracted {} signatures", featuresWithSignatureCount);
    }

    private static void writeFusion(@NotNull String fusionTsv, @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry)
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(fusionTsv));
        writer.write("Fusion" + "\t" + "Fusion_event" + "\n");

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {

            for (Map.Entry<Feature, FusionAnnotation> featureResult : entry.getValue().fusionsPerFeature().entrySet()) {
                FusionAnnotation geneFusionForFeature = featureResult.getValue();

                writer.write(
                        geneFusionForFeature.fusion() + "\t" + geneFusionForFeature.fusionEvent() + "\n");
            }
        }
        writer.close();
    }

    private static void writeRanges(@NotNull String rangesTsv, @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry)
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(rangesTsv));
        writer.write("Gene" + "\t" + "chromosome" + "\t" + "start" + "\t" + "end" + "\t" + "event" + "\n");

        for (Map.Entry<ViccEntry, ViccExtractionResult> entry : resultsPerEntry.entrySet()) {
            for (Map.Entry<Feature, GeneRangeAnnotation> featureResult : entry.getValue().geneRangesPerFeature().entrySet()) {
                GeneRangeAnnotation geneRangeForFeature = featureResult.getValue();

                writer.write(
                        geneRangeForFeature.gene() + "\t" + geneRangeForFeature.chromosome() + "\t" + geneRangeForFeature.start() + "\t"
                                + geneRangeForFeature.end() + "\t" + geneRangeForFeature.event() + "\n");
            }
        }
        writer.close();
    }

    private static void writeHotspots(@NotNull String hotspotVcf, @NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(hotspotVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        writer.writeHeader(header);

        for (Map.Entry<VariantHotspot, HotspotAnnotation> entry : convertAndSort(resultsPerEntry).entrySet()) {
            VariantHotspot hotspot = entry.getKey();
            HotspotAnnotation annotation = entry.getValue();
            List<Allele> hotspotAlleles = buildAlleles(hotspot);

            VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                    .source("VICC")
                    .chr(hotspot.chromosome())
                    .start(hotspot.position())
                    .alleles(hotspotAlleles)
                    .computeEndFromAlleles(hotspotAlleles, (int) hotspot.position())
                    .attribute("sources", buildSourcesString(annotation.sources()))
                    .attribute("feature",
                            ProteinKeyFormatter.toProteinKey(annotation.gene(), annotation.transcript(), annotation.proteinAnnotation()))
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

    @VisibleForTesting
    @NotNull
    static String buildSourcesString(@NotNull Set<String> sources) {
        StringJoiner sourceJoiner = new StringJoiner(",");
        for (String source : sources) {
            sourceJoiner.add(source);
        }
        return sourceJoiner.toString();
    }

    @NotNull
    private static Map<VariantHotspot, HotspotAnnotation> convertAndSort(@NotNull Map<ViccEntry, ViccExtractionResult> resultsPerEntry) {
        Map<VariantHotspot, HotspotAnnotation> convertedMap = Maps.newTreeMap(new VariantHotspotComparator());
        for (Map.Entry<ViccEntry, ViccExtractionResult> entryResult : resultsPerEntry.entrySet()) {
            ViccEntry entry = entryResult.getKey();
            for (Map.Entry<Feature, List<VariantHotspot>> featureResult : entryResult.getValue().hotspotsPerFeature().entrySet()) {
                Feature feature = featureResult.getKey();
                for (VariantHotspot hotspot : featureResult.getValue()) {
                    HotspotAnnotation annotation = convertedMap.get(hotspot);
                    if (annotation != null) {
                        if (annotation.sources().contains(entry.source().display())) {
                            checkForDuplicateHotspotOnDifferentProteinAnnotation(entry, annotation, feature);

                            // We try to override transcript in case we got a non-null one.
                            String bestTranscript = annotation.transcript() == null ? entry.transcriptId() : annotation.transcript();
                            annotation = new HotspotAnnotation(annotation.sources(),
                                    annotation.gene(),
                                    bestTranscript,
                                    annotation.proteinAnnotation());

                        } else {
                            Set<String> newSources = Sets.newHashSet(annotation.sources());
                            newSources.add(entry.source().display());
                            String bestTranscript = annotation.transcript() == null ? entry.transcriptId() : annotation.transcript();
                            annotation =
                                    new HotspotAnnotation(newSources, annotation.gene(), bestTranscript, annotation.proteinAnnotation());
                        }
                    } else {
                        annotation = new HotspotAnnotation(Sets.newHashSet(entry.source().display()),
                                feature.geneSymbol(),
                                entry.transcriptId(),
                                feature.proteinAnnotation());
                    }

                    convertedMap.put(hotspot, annotation);
                }
            }
        }

        return convertedMap;
    }

    private static void checkForDuplicateHotspotOnDifferentProteinAnnotation(@NotNull ViccEntry entry,
            @NotNull HotspotAnnotation annotation, @NotNull Feature feature) {
        String transcript = entry.transcriptId();
        if (annotation.gene().equals(feature.geneSymbol()) && transcript != null && transcript.equals(annotation.transcript())
                && !feature.proteinAnnotation().equals(annotation.proteinAnnotation())) {
            String existingKey =
                    ProteinKeyFormatter.toProteinKey(annotation.gene(), annotation.transcript(), annotation.proteinAnnotation());
            String newKey = ProteinKeyFormatter.toProteinKey(feature.geneSymbol(), transcript, feature.proteinAnnotation());
//                        LOGGER.warn("Hotspot already exists for '{}' under a different annotation. " + "Existing key = '{}'. New key = '{}'",
//                                entry.source().display(),
//                                existingKey,
//                                newKey);
        }
    }
}
