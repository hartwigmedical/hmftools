package com.hartwig.hmftools.serve;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.serve.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.sources.ExtractionOutput;
import com.hartwig.hmftools.serve.sources.ImmutableExtractionOutput;
import com.hartwig.hmftools.serve.sources.docm.DocmEntry;
import com.hartwig.hmftools.serve.sources.docm.DocmExtractor;
import com.hartwig.hmftools.serve.sources.docm.DocmReader;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigEntry;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigExtractor;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigFileReader;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionExtractor;
import com.hartwig.hmftools.serve.sources.iclusion.IclusionReader;
import com.hartwig.hmftools.serve.sources.vicc.ViccExtractor;
import com.hartwig.hmftools.serve.sources.vicc.ViccExtractorFactory;
import com.hartwig.hmftools.serve.sources.vicc.ViccReader;
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

        ProteinResolver proteinResolver = buildProteinResolver(config);
        List<DriverGene> driverGenes = readDriverGenes(config);

        List<ExtractionOutput> extractions = Lists.newArrayList();
        extractions.add(extractViccKnowledge(config.viccJson(), proteinResolver, driverGenes));
        extractions.add(extractIclusionKnowledge(config.iClusionTrialTsv()));
        extractions.add(extractDocmKnowledge(config.docmTsv(), proteinResolver));
        extractions.add(extractHartwigCohortKnowledge(config.hartwigCohortTsv(), proteinResolver, !config.skipHotspotResolving()));
        extractions.add(extractHartwigCuratedKnowledge(config.hartwigCuratedTsv(), proteinResolver, !config.skipHotspotResolving()));

        LOGGER.info("Done with {}", extractions);
    }

    @NotNull
    private static ProteinResolver buildProteinResolver(@NotNull ServeConfig config) throws FileNotFoundException {
        if (config.skipHotspotResolving()) {
            LOGGER.info("Creating dummy protein resolver");
            return ProteinResolverFactory.dummy();
        } else {
            LOGGER.info("Creating transvar protein resolver with ref genome version '{}' using FASTA {}",
                    config.refGenomeVersion(),
                    config.refGenomeFastaFile());
            return ProteinResolverFactory.transvarWithRefGenome(config.refGenomeVersion(), config.refGenomeFastaFile());
        }
    }

    @NotNull
    private static List<DriverGene> readDriverGenes(@NotNull ServeConfig config) throws IOException {
        LOGGER.info("Reading driver genes from {}", config.driverGeneTsv());
        List<DriverGene> driverGenes = DriverGeneFile.read(config.driverGeneTsv());
        LOGGER.info(" Read {} driver gene entries", driverGenes.size());
        return driverGenes;
    }

    @NotNull
    private static ExtractionOutput extractViccKnowledge(@NotNull String viccJson, @NotNull ProteinResolver proteinResolver,
            @NotNull List<DriverGene> driverGenes) throws IOException {
        List<ViccEntry> entries = ViccReader.readAndCurateRelevantEntries(viccJson, VICC_SOURCES_TO_INCLUDE, null);

        ViccExtractor extractor = ViccExtractorFactory.buildViccExtractor(proteinResolver, driverGenes);
        return extractor.extractFromViccEntries(entries);
    }

    @NotNull
    private static ExtractionOutput extractIclusionKnowledge(@NotNull String iClusionTrialTsv) throws IOException {
        List<IclusionTrial> trials = IclusionReader.readAndCurate(iClusionTrialTsv);

        IclusionExtractor extractor = new IclusionExtractor();
        return extractor.extractFromIclusionTrials(trials);
    }

    @NotNull
    private static ExtractionOutput extractDocmKnowledge(@NotNull String docmTsv, @NotNull ProteinResolver proteinResolver)
            throws IOException {
        List<DocmEntry> entries = DocmReader.readAndCurate(docmTsv);

        DocmExtractor extractor = new DocmExtractor(proteinResolver);
        List<KnownHotspot> hotspots = extractor.extractFromDocmEntries(entries);
        return ImmutableExtractionOutput.builder().knownHotspots(hotspots).build();
    }

    @NotNull
    private static ExtractionOutput extractHartwigCohortKnowledge(@NotNull String hartwigCohortTsv,
            @NotNull ProteinResolver proteinResolver, boolean addExplicitHotspots) throws IOException {
        LOGGER.info("Reading Hartwig Cohort TSV from '{}'", hartwigCohortTsv);
        List<HartwigEntry> entries = HartwigFileReader.read(hartwigCohortTsv);
        LOGGER.info(" Read {} entries", entries.size());

        HartwigExtractor extractor = new HartwigExtractor(Knowledgebase.HARTWIG_COHORT, proteinResolver, addExplicitHotspots);

        List<KnownHotspot> hotspots =  extractor.extractFromHartwigEntries(entries);
        return ImmutableExtractionOutput.builder().knownHotspots(hotspots).build();
    }

    @NotNull
    private static ExtractionOutput extractHartwigCuratedKnowledge(@NotNull String hartwigCuratedTsv,
            @NotNull ProteinResolver proteinResolver, boolean addExplicitHotspots) throws IOException {
        LOGGER.info("Reading Hartwig Curated TSV from '{}'", hartwigCuratedTsv);
        List<HartwigEntry> entries = HartwigFileReader.read(hartwigCuratedTsv);
        LOGGER.info(" Read {} entries", entries.size());

        HartwigExtractor extractor = new HartwigExtractor(Knowledgebase.HARTWIG_CURATED, proteinResolver, addExplicitHotspots);

        List<KnownHotspot> hotspots = extractor.extractFromHartwigEntries(entries);
        return ImmutableExtractionOutput.builder().knownHotspots(hotspots).build();
    }
}
