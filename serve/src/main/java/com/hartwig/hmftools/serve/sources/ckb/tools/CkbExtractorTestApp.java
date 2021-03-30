package com.hartwig.hmftools.serve.sources.ckb.tools;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import javax.lang.model.element.ElementVisitor;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.JsonDatabaseToCkbEntryConverter;
import com.hartwig.hmftools.ckb.classification.CkbClassificationConfig;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.CkbJsonReader;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.curation.DoidLookupFactory;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.refgenome.ImmutableRefGenomeResource;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.ckb.CkbExtractor;
import com.hartwig.hmftools.serve.sources.ckb.CkbExtractorFactory;
import com.hartwig.hmftools.serve.sources.ckb.CkbReader;
import com.hartwig.hmftools.serve.sources.ckb.CkbUtils;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractor;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractorFactory;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionUtil;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;

public class CkbExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractorTestApp.class);
    private static final String FIELD_DELIMITER = "\t";

    public static void main(String[] args) throws IOException {
        String ckbDir = "/data/common/dbs/ckb/210326_flex_dump";
        String outputDir = System.getProperty("user.home") + "/tmp";

        CkbJsonDatabase ckbJsonDatabase = CkbJsonReader.read(ckbDir);
        List<CkbEntry> allCkbEntries = JsonDatabaseToCkbEntryConverter.convert(ckbJsonDatabase);
        List<CkbEntry> filteredAndcurateCkbEntries = CkbReader.filterAndCurateRelevantEntries(allCkbEntries);

        EventClassifierConfig config = CkbClassificationConfig.build();
        EventClassifier classifier = EventClassifierFactory.buildClassifier(config);

        // TODO 1. Make sure every entry has correct event type.

        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("gene").add("event").add("type").toString();
        lines.add(header);

        EventType type = EventType.UNKNOWN;

        for (CkbEntry entry : filteredAndcurateCkbEntries) {
            String gene = entry.variants().get(0).gene().geneSymbol();

            String profileName;
            if (entry.variants().size() > 1) {
                type = EventType.COMBINED;
                profileName = entry.profileName();
            } else {
                if (entry.variants().get(0).variant().equals("fusion") && entry.variants().get(0).impact() != null && entry.variants()
                        .get(0)
                        .impact()
                        .equals("fusion")) {
                    profileName = "fusion promisuous";
                } else if (entry.variants().get(0).impact() != null && entry.variants().get(0).impact().equals("fusion")) {
                    profileName = entry.variants().get(0).variant().replaceAll("\\s+","") + " fusion";
                } else if (entry.variants().get(0).variant().contains("exon")) {
                    profileName = entry.variants().get(0).variant().replace("exon", "exon ");
                }
                else {
                    profileName = entry.variants().get(0).variant() ;
                }

                type = classifier.determineType(gene, profileName);
            }

            lines.add(new StringJoiner(FIELD_DELIMITER).add(gene).add(profileName).add(type.toString()).toString());
        }
        Files.write(new File("/data/common/dbs/serve/pilot_output/events.tsv").toPath(), lines);

        // TODO 2. Make sure every event is extracted correctly using EventExtractor
        String driverGeneTsvPath = "/data/common/dbs/driver_gene_panel/DriverGenePanel.hg38.tsv";
        String knownFusionFilePath = "/data/common/dbs/fusions/known_fusion_data.csv";
        String missingDoidMappingTsv = "/data/common/dbs/serve/curation/missing_doids_mapping.tsv";
        String fastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh38.no.alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";
        RefGenomeVersion refGenomeVersion = RefGenomeVersion.V38;
        Map<String, HmfTranscriptRegion> allGenesMap = HmfGenePanelSupplier.allGenesMap38();
        ProteinResolver proteinResolver = ProteinResolverFactory.dummy();


        DoidLookup doidLookup = DoidLookupFactory.buildFromConfigTsv(missingDoidMappingTsv);
        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsvPath);
        LOGGER.debug(" Read {} driver genes from {}", driverGenes.size(), driverGeneTsvPath);

        KnownFusionCache fusionCache = new KnownFusionCache();
        if (!fusionCache.loadFile(knownFusionFilePath)) {
            throw new IllegalStateException("Could not load known fusion cache from " + knownFusionFilePath);
        }
        LOGGER.debug(" Read {} known fusions from {}", fusionCache.getData().size(), knownFusionFilePath);
        RefGenomeResource refGenomeResource = ImmutableRefGenomeResource.builder()
                .fastaFile(fastaFile)
                .driverGenes(driverGenes)
                .knownFusionCache(fusionCache)
                .canonicalTranscriptPerGeneMap(allGenesMap)
                .proteinResolver(proteinResolver)
                .build();

        CkbExtractor extractor = CkbExtractorFactory.buildCkbExtractor(config, refGenomeResource, doidLookup);

        // TODO 3. Create ActionableEvents for all relevant entries.
        ExtractionResult result = extractor.extract(filteredAndcurateCkbEntries, type);

        CkbUtils.printExtractionResults(result, type);


        // IclusionUtil.printIclusionResult(result);
     //   IclusionUtil.writeIclusionMutationTypes(iclusionMutationTsv, trials);

        new ExtractionResultWriter(outputDir, refGenomeVersion).write(result);

    }
}
