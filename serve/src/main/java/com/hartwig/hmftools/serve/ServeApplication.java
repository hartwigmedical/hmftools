package com.hartwig.hmftools.serve;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.iclusion.classification.IclusionClassificationConfig;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.serve.extraction.ExtractionFunctions;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ExtractionResultWriter;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.sources.docm.DocmEntry;
import com.hartwig.hmftools.serve.sources.docm.DocmExtractor;
import com.hartwig.hmftools.serve.sources.docm.DocmReader;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigEntry;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigExtractor;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigFileReader;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractor;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractorFactory;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionReader;
import com.hartwig.hmftools.serve.sources.vicc.ViccExtractor;
import com.hartwig.hmftools.serve.sources.vicc.ViccExtractorFactory;
import com.hartwig.hmftools.serve.sources.vicc.ViccReader;
import com.hartwig.hmftools.vicc.annotation.ViccClassificationConfig;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ServeApplication {

    private static final Logger LOGGER = LogManager.getLogger(ServeApplication.class);

    private static final Set<ViccSource> VICC_SOURCES_TO_INCLUDE =
            Sets.newHashSet(ViccSource.CIVIC, ViccSource.JAX, ViccSource.ONCOKB, ViccSource.CGI);

    private static final String VERSION = ServeApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE v{}", VERSION);

        Options options = ServeConfig.createOptions();

        ServeConfig config = null;
        try {
            config = ServeConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("SERVE", options);
            System.exit(1);
        }

        List<DriverGene> driverGenes = readDriverGenesFromFile(config.driverGeneTsv());
        Map<String, HmfTranscriptRegion> allGenesMap = selectAllGeneMap(config);
        ProteinResolver proteinResolver = buildProteinResolver(config, allGenesMap);

        List<ExtractionResult> extractions = Lists.newArrayList();
        extractions.add(extractViccKnowledge(config.viccJson(), proteinResolver, driverGenes, allGenesMap));
        extractions.add(extractIclusionKnowledge(config.iClusionTrialTsv(), proteinResolver, driverGenes, allGenesMap));
        extractions.add(extractDocmKnowledge(config.docmTsv(), proteinResolver));
        extractions.add(extractHartwigCohortKnowledge(config.hartwigCohortTsv(), proteinResolver, !config.skipHotspotResolving()));
        extractions.add(extractHartwigCuratedKnowledge(config.hartwigCuratedTsv(), proteinResolver, !config.skipHotspotResolving()));

        evaluateProteinResolver(proteinResolver);

        ExtractionResultWriter writer = new ExtractionResultWriter(config.outputDir(), config.refGenomeVersion());
        writer.write(ExtractionFunctions.merge(extractions));

        LOGGER.info("Done!");
    }

    @NotNull
    private static ProteinResolver buildProteinResolver(@NotNull ServeConfig config,
            @NotNull Map<String, HmfTranscriptRegion> transcriptsPerGeneMap) throws FileNotFoundException {
        if (config.skipHotspotResolving()) {
            LOGGER.info("Creating dummy protein resolver");
            return ProteinResolverFactory.dummy();
        } else {
            LOGGER.info("Creating transvar protein resolver with ref genome version '{}' using FASTA {}",
                    config.refGenomeVersion(),
                    config.refGenomeFastaFile());
            return ProteinResolverFactory.transvarWithRefGenome(config.refGenomeVersion(),
                    config.refGenomeFastaFile(),
                    transcriptsPerGeneMap);
        }
    }

    private static void evaluateProteinResolver(@NotNull ProteinResolver proteinResolver) {
        Set<String> unresolvedProteinAnnotations = proteinResolver.unresolvedProteinAnnotations();
        if (!unresolvedProteinAnnotations.isEmpty()) {
            LOGGER.warn("Protein resolver could not resolve {} protein annotations", unresolvedProteinAnnotations.size());
            for (String unresolvedProteinAnnotation : unresolvedProteinAnnotations) {
                LOGGER.warn("Protein resolver could not resolve protein annotation '{}'", unresolvedProteinAnnotation);
            }
        } else {
            LOGGER.debug("Protein resolver observed no issues when resolving hotspots");
        }
    }

    @NotNull
    private static List<DriverGene> readDriverGenesFromFile(@NotNull String driverGeneTsv) throws IOException {
        LOGGER.info("Reading driver genes from {}", driverGeneTsv);
        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsv);
        LOGGER.info(" Read {} driver gene entries", driverGenes.size());
        return driverGenes;
    }

    @NotNull
    private static Map<String, HmfTranscriptRegion> selectAllGeneMap(@NotNull ServeConfig config) {
        if (config.refGenomeVersion() == RefGenomeVersion.V37) {
            return HmfGenePanelSupplier.allGenesMap37();
        }

        LOGGER.warn("Ref genome version not supported: '{}'. Reverting to HG19.", config.refGenomeVersion());
        return HmfGenePanelSupplier.allGenesMap37();
    }

    @NotNull
    private static ExtractionResult extractViccKnowledge(@NotNull String viccJson, @NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes, @NotNull Map<String, HmfTranscriptRegion> allGenesMap) throws IOException {
        List<ViccEntry> entries = ViccReader.readAndCurateRelevantEntries(viccJson, VICC_SOURCES_TO_INCLUDE, null);

        EventClassifierConfig config = ViccClassificationConfig.build();
        ViccExtractor extractor = ViccExtractorFactory.buildViccExtractor(config, proteinResolver, driverGenes, allGenesMap);
        LOGGER.info("Running VICC knowledge extraction");
        return extractor.extract(entries);
    }

    @NotNull
    private static ExtractionResult extractIclusionKnowledge(@NotNull String iClusionTrialTsv, @NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes, @NotNull Map<String, HmfTranscriptRegion> allGenesMap) throws IOException {
        List<IclusionTrial> trials = IclusionReader.readAndCurate(iClusionTrialTsv);

        EventClassifierConfig config = IclusionClassificationConfig.build();
        IclusionExtractor extractor = IclusionExtractorFactory.buildIclusionExtractor(config, proteinResolver, driverGenes, allGenesMap);
        LOGGER.info("Running iClusion knowledge extraction");
        return extractor.extract(trials);
    }

    @NotNull
    private static ExtractionResult extractDocmKnowledge(@NotNull String docmTsv, @NotNull ProteinResolver proteinResolver)
            throws IOException {
        List<DocmEntry> entries = DocmReader.readAndCurate(docmTsv);

        DocmExtractor extractor = new DocmExtractor(proteinResolver);
        LOGGER.info("Running DoCM knowledge extraction");
        return extractor.extract(entries);
    }

    @NotNull
    private static ExtractionResult extractHartwigCohortKnowledge(@NotNull String hartwigCohortTsv,
            @NotNull ProteinResolver proteinResolver, boolean addExplicitHotspots) throws IOException {
        LOGGER.info("Reading Hartwig Cohort TSV from '{}'", hartwigCohortTsv);
        List<HartwigEntry> entries = HartwigFileReader.read(hartwigCohortTsv);
        LOGGER.info(" Read {} entries", entries.size());

        HartwigExtractor extractor = new HartwigExtractor(Knowledgebase.HARTWIG_COHORT, proteinResolver, addExplicitHotspots);
        LOGGER.info("Running Hartwig Cohort knowledge extraction");
        return extractor.extract(entries);
    }

    @NotNull
    private static ExtractionResult extractHartwigCuratedKnowledge(@NotNull String hartwigCuratedTsv,
            @NotNull ProteinResolver proteinResolver, boolean addExplicitHotspots) throws IOException {
        LOGGER.info("Reading Hartwig Curated TSV from '{}'", hartwigCuratedTsv);
        List<HartwigEntry> entries = HartwigFileReader.read(hartwigCuratedTsv);
        LOGGER.info(" Read {} entries", entries.size());

        HartwigExtractor extractor = new HartwigExtractor(Knowledgebase.HARTWIG_CURATED, proteinResolver, addExplicitHotspots);
        LOGGER.info("Running Hartwig Curated knowledge extraction");
        return extractor.extract(entries);
    }
}
