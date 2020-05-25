package com.hartwig.hmftools.serve.vicc;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.copynumber.KnownAmplificationDeletion;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class ViccExtractorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccExtractorTestApplication.class);

    private static final boolean RUN_ON_SERVER = false;

    public static void main(String[] args) throws IOException, InterruptedException {
        Configurator.setRootLevel(Level.DEBUG);

        String viccJsonPath;
        String refGenomeFastaFile;

        if (RUN_ON_SERVER) {
            viccJsonPath = "/data/common/dbs/vicc/all.json";
            refGenomeFastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
        } else {
            viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";
            refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
        }

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;

        String source = "oncokb";
        LOGGER.info("Reading VICC json from {} with source '{}'", viccJsonPath, source);
        List<ViccEntry> viccEntries = ViccJsonReader.readSingleKnowledgebase(viccJsonPath, source);
        LOGGER.info(" Read {} entries", viccEntries.size());

        HotspotExtractor hotspotExtractor = HotspotExtractor.withRefGenome(refGenomeVersion, refGenomeFastaFile);
        CopyNumberExtractor copyNumberExtractor = new CopyNumberExtractor();
        FusionExtractor fusionExtractor = new FusionExtractor();

        Map<Feature, List<VariantHotspot>> allHotspotsPerFeature = Maps.newHashMap();
        Map<Feature, KnownAmplificationDeletion> allKnownAmpsDelsPerFeature = Maps.newHashMap();
        Map<Feature, String> allKnownFusions = Maps.newHashMap();
        for (ViccEntry viccEntry : viccEntries) {
            allHotspotsPerFeature.putAll(hotspotExtractor.extractHotspots(viccEntry));
            allKnownAmpsDelsPerFeature.putAll(copyNumberExtractor.extractKnownAmplificationsDeletions(viccEntry));
            allKnownFusions.putAll(fusionExtractor.extractKnownFusions(viccEntry));
        }

        List<Feature> featuresWithoutGenomicEvents = Lists.newArrayList();
        int totalFeatureCount = 0;
        for (ViccEntry viccEntry : viccEntries) {
            for (Feature feature : viccEntry.features()) {
                if (!allHotspotsPerFeature.containsKey(feature) && !allKnownAmpsDelsPerFeature.containsKey(feature)
                        && !allKnownFusions.containsKey(feature)) {
                    featuresWithoutGenomicEvents.add(feature);
                }
                totalFeatureCount++;
            }
        }
        LOGGER.info("Done extraction from '{}'", source);
        LOGGER.info(" Extraction performed on {} features from {} entries", totalFeatureCount, viccEntries.size());
        LOGGER.info(" Extracted {} hotspots for {} features", valuesCount(allHotspotsPerFeature), allHotspotsPerFeature.size());
        LOGGER.info(" Extracted {} known amps and dels", allKnownAmpsDelsPerFeature.size());
        LOGGER.info(" Extracted {} fusions", allKnownFusions.size());

        LOGGER.info("Could not resolve hotspots for {} features", hotspotExtractor.unresolvableFeatures().size());
        for (String feature : hotspotExtractor.unresolvableFeatures()) {
            LOGGER.debug(" {}", feature);
        }

        LOGGER.info("No genomic events found for {} features", featuresWithoutGenomicEvents.size());
        for (Feature feature : featuresWithoutGenomicEvents) {
            LOGGER.debug(" {}", feature);
        }
    }

    private static <T, Y> int valuesCount(@NotNull Map<T, List<Y>> map) {
        int valuesCount = 0;
        for (Map.Entry<T, List<Y>> entry : map.entrySet()) {
            valuesCount += entry.getValue().size();
        }
        return valuesCount;
    }
}
