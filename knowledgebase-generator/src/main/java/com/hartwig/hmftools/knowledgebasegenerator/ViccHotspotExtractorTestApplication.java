package com.hartwig.hmftools.knowledgebasegenerator;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.knowledgebasegenerator.hotspot.HotspotExtractor;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ViccHotspotExtractorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccHotspotExtractorTestApplication.class);

    public static void main(String[] args) throws IOException, InterruptedException {
        String viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";
        String refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;

        String source = "oncokb";
        LOGGER.info("Reading VICC json from {} with source '{}'", viccJsonPath, source);
        List<ViccEntry> viccEntries = ViccJsonReader.readSingleKnowledgebase(viccJsonPath, source);
        LOGGER.info("Read {} entries", viccEntries.size());

        HotspotExtractor hotspotExtractor = HotspotExtractor.withRefGenome(refGenomeVersion, refGenomeFastaFile);
        for (ViccEntry viccEntry : viccEntries) {
            hotspotExtractor.extractHotspots(viccEntry);
        }
    }
}
