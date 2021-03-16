package com.hartwig.hmftools.serve;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.iclusion.classification.IclusionClassificationConfig;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.ExtractionFunctions;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ServeAlgo {

    private static final Logger LOGGER = LogManager.getLogger(ServeAlgo.class);

    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;
    @NotNull
    private final Map<String, HmfTranscriptRegion> allGenesMap;
    @NotNull
    private final ProteinResolver proteinResolver;
    @NotNull
    private final DoidLookup missingDoidLookup;

    public ServeAlgo(@NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache,
            @NotNull final Map<String, HmfTranscriptRegion> allGenesMap, @NotNull final ProteinResolver proteinResolver,
            @NotNull final DoidLookup missingDoidLookup) {
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
        this.allGenesMap = allGenesMap;
        this.proteinResolver = proteinResolver;
        this.missingDoidLookup = missingDoidLookup;
    }

    @NotNull
    public ExtractionResult run(@NotNull ServeConfig config)
            throws IOException {
        List<ExtractionResult> extractions = Lists.newArrayList();
        if (config.useVicc()) {
            extractions.add(extractViccKnowledge(config.viccJson(), config.viccSources()));
        }

        if (config.useIclusion()) {
            extractions.add(extractIclusionKnowledge(config.iClusionTrialTsv()));
        }

        if (config.useDocm()) {
            extractions.add(extractDocmKnowledge(config.docmTsv()));
        }

        if (config.useHartwigCohort()) {
            extractions.add(extractHartwigCohortKnowledge(config.hartwigCohortTsv(), !config.skipHotspotResolving()));
        }

        if (config.useHartwigCurated()) {
            extractions.add(extractHartwigCuratedKnowledge(config.hartwigCuratedTsv(), !config.skipHotspotResolving()));
        }

        evaluateProteinResolver(proteinResolver);
        missingDoidLookup.reportUnusedMappings();

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private ExtractionResult extractViccKnowledge(@NotNull String viccJson, @NotNull Set<ViccSource> viccSources) throws IOException {
        List<ViccEntry> entries = ViccReader.readAndCurateRelevantEntries(viccJson, viccSources, null);

        EventClassifierConfig config = ViccClassificationConfig.build();
        ViccExtractor extractor = ViccExtractorFactory.buildViccExtractor(config,
                proteinResolver,
                driverGenes,
                knownFusionCache,
                allGenesMap,
                missingDoidLookup);

        LOGGER.info("Running VICC knowledge extraction");
        return extractor.extract(entries);
    }

    @NotNull
    private ExtractionResult extractIclusionKnowledge(@NotNull String iClusionTrialTsv) throws IOException {
        List<IclusionTrial> trials = IclusionReader.readAndCurate(iClusionTrialTsv);

        EventClassifierConfig config = IclusionClassificationConfig.build();
        IclusionExtractor extractor = IclusionExtractorFactory.buildIclusionExtractor(config,
                proteinResolver,
                driverGenes,
                knownFusionCache,
                allGenesMap,
                missingDoidLookup);

        LOGGER.info("Running iClusion knowledge extraction");
        return extractor.extract(trials);
    }

    @NotNull
    private ExtractionResult extractDocmKnowledge(@NotNull String docmTsv) throws IOException {
        List<DocmEntry> entries = DocmReader.readAndCurate(docmTsv);

        DocmExtractor extractor = new DocmExtractor(proteinResolver);
        LOGGER.info("Running DoCM knowledge extraction");
        return extractor.extract(entries);
    }

    @NotNull
    private ExtractionResult extractHartwigCohortKnowledge(@NotNull String hartwigCohortTsv, boolean addExplicitHotspots)
            throws IOException {
        LOGGER.info("Reading Hartwig Cohort TSV from '{}'", hartwigCohortTsv);
        List<HartwigEntry> entries = HartwigFileReader.read(hartwigCohortTsv);
        LOGGER.info(" Read {} entries", entries.size());

        HartwigExtractor extractor = new HartwigExtractor(Knowledgebase.HARTWIG_COHORT, proteinResolver, addExplicitHotspots);
        LOGGER.info("Running Hartwig Cohort knowledge extraction");
        return extractor.extract(entries);
    }

    @NotNull
    private ExtractionResult extractHartwigCuratedKnowledge(@NotNull String hartwigCuratedTsv, boolean addExplicitHotspots)
            throws IOException {
        LOGGER.info("Reading Hartwig Curated TSV from '{}'", hartwigCuratedTsv);
        List<HartwigEntry> entries = HartwigFileReader.read(hartwigCuratedTsv);
        LOGGER.info(" Read {} entries", entries.size());

        HartwigExtractor extractor = new HartwigExtractor(Knowledgebase.HARTWIG_CURATED, proteinResolver, addExplicitHotspots);
        LOGGER.info("Running Hartwig Curated knowledge extraction");
        return extractor.extract(entries);
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
}
