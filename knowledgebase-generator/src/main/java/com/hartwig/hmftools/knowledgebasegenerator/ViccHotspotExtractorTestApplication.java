package com.hartwig.hmftools.knowledgebasegenerator;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.knowledgebasegenerator.hotspot.HotspotExtractor;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class ViccHotspotExtractorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccHotspotExtractorTestApplication.class);

    private static final boolean RUN_ON_SERVER = true;

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
        LOGGER.info("Read {} entries", viccEntries.size());

        HotspotExtractor hotspotExtractor = HotspotExtractor.withRefGenome(refGenomeVersion, refGenomeFastaFile);
        Map<String, List<VariantHotspot>> allHotspotsPerFeature = Maps.newHashMap();
        for (ViccEntry viccEntry : viccEntries) {
            mergeIntoExistingMap(allHotspotsPerFeature, hotspotExtractor.extractHotspots(viccEntry));
        }

        LOGGER.info("Done extracting {} hotspots for {} features", valuesCount(allHotspotsPerFeature), allHotspotsPerFeature.size());

        LOGGER.info("Printing unresolvable features");
        for (String feature : hotspotExtractor.unresolvableFeatures()) {
            LOGGER.info(feature);
        }
    }

    @VisibleForTesting
    static <T, Y> void mergeIntoExistingMap(@NotNull Map<T, List<Y>> map, @NotNull Map<T, List<Y>> mapToMergeIn) {
        for (Map.Entry<T, List<Y>> entry : mapToMergeIn.entrySet()) {
            if (map.containsKey(entry.getKey())) {
                List<Y> currentEntries = map.get(entry.getKey());
                currentEntries.addAll(entry.getValue());
                map.put(entry.getKey(), currentEntries);
            } else {
                map.put(entry.getKey(), entry.getValue());
            }
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
