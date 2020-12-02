package com.hartwig.hmftools.serve.hotspot.tools;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.hotspot.HotspotVCF;
import com.hartwig.hmftools.serve.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.sources.docm.DocmEntry;
import com.hartwig.hmftools.serve.sources.docm.DocmExtractor;
import com.hartwig.hmftools.serve.sources.docm.DocmReader;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigEntry;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigExtractor;
import com.hartwig.hmftools.serve.sources.hartwig.HartwigFileReader;
import com.hartwig.hmftools.serve.sources.vicc.ViccExtractor;
import com.hartwig.hmftools.serve.sources.vicc.ViccExtractorFactory;
import com.hartwig.hmftools.serve.sources.vicc.ViccReader;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class ServeHotspotGenerator {

    private static final Logger LOGGER = LogManager.getLogger(ServeHotspotGenerator.class);

    private static final Set<ViccSource> VICC_SOURCES_TO_INCLUDE = Sets.newHashSet(ViccSource.CIVIC, ViccSource.CGI);
    private static final Integer MAX_VICC_ENTRIES = null;

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String refGenomeFastaFile;
        String serveSourceDir;

        boolean generateHotspots;
        String hotspotVcf = null;

        if (hostname.toLowerCase().contains("datastore")) {
            refGenomeFastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            serveSourceDir = "/data/common/dbs/serve/static_sources";

            generateHotspots = true;
            hotspotVcf = System.getProperty("user.home") + "/tmp/hotspotsServe.vcf";
        } else {
            refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            serveSourceDir = System.getProperty("user.home") + "/hmf/projects/serve";

            generateHotspots = false;
        }

        String viccJson = serveSourceDir + "/vicc/all.json";
        String docmTsv = serveSourceDir + "/docm/docm_v3.2.tsv";
        String hartwigCohortTsv = serveSourceDir + "/hartwig/hartwig_cohort.tsv";
        String hartwigCuratedTsv = serveSourceDir + "/hartwig/hartwig_curated.tsv";

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        LOGGER.debug("Configured '{}' as the reference fasta path", refGenomeFastaFile);
        LOGGER.debug("Configured '{}' as the VICC json path", viccJson);
        LOGGER.debug("Configured '{}' as the DoCM TSV path", docmTsv);
        LOGGER.debug("Configured '{}' as the Hartwig Cohort TSV path", hartwigCohortTsv);
        LOGGER.debug("Configured '{}' as the Hartwig Curated TSV path", hartwigCuratedTsv);
        LOGGER.debug("Configured '{}' as the hotspot output VCF", hotspotVcf);
        LOGGER.debug("Configured '{}' for generating hotspots yes/no", generateHotspots);

        ProteinResolver proteinResolver = generateHotspots
                ? ProteinResolverFactory.transvarWithRefGenome(refGenomeVersion, refGenomeFastaFile, HmfGenePanelSupplier.allGenesMap37())
                : ProteinResolverFactory.dummy();

        List<KnownHotspot> hartwigCohortHotspots = hartwigCohortHotspots(hartwigCohortTsv, proteinResolver, generateHotspots);
        List<KnownHotspot> hartwigCuratedHotspots = hartwigCuratedHotspots(hartwigCuratedTsv, proteinResolver, generateHotspots);
        List<KnownHotspot> viccHotspots = viccHotspots(viccJson, proteinResolver);
        List<KnownHotspot> docmHotspots = docmHotspots(docmTsv, proteinResolver);

        LOGGER.info("Merging {} VICC hotspots with {} DoCM hotspots and {} Hartwig Cohort hotspots and {} Hartwig Curated hotspots",
                viccHotspots.size(),
                docmHotspots.size(),
                hartwigCohortHotspots.size(),
                hartwigCuratedHotspots.size());
        List<KnownHotspot> allHotspots = Lists.newArrayList();
        allHotspots.addAll(viccHotspots);
        allHotspots.addAll(docmHotspots);
        allHotspots.addAll(hartwigCohortHotspots);
        allHotspots.addAll(hartwigCuratedHotspots);

        List<KnownHotspot> mergedConsolidatedHotspots = HotspotFunctions.consolidate(allHotspots);

        if (generateHotspots && hotspotVcf != null) {
            HotspotVCF.write(hotspotVcf, mergedConsolidatedHotspots);

            Set<String> unresolvedProteinAnnotations = proteinResolver.unresolvedProteinAnnotations();
            if (!unresolvedProteinAnnotations.isEmpty()) {
                LOGGER.warn("Protein resolver could not resolve {} protein annotations", unresolvedProteinAnnotations.size());
                for (String unresolvedProteinAnnotation : unresolvedProteinAnnotations) {
                    LOGGER.warn("Protein resolver could not resolve protein annotation '{}'", unresolvedProteinAnnotation);
                }
            } else {
                LOGGER.info("Protein resolver could resolve hotspots for every protein annotation");
            }
        }
    }

    @NotNull
    private static List<KnownHotspot> viccHotspots(@NotNull String viccJson, @NotNull ProteinResolver proteinResolver) throws IOException {
        List<ViccEntry> viccEntries = ViccReader.readAndCurateRelevantEntries(viccJson, VICC_SOURCES_TO_INCLUDE, MAX_VICC_ENTRIES);

        ViccExtractor viccExtractor =
                ViccExtractorFactory.buildViccExtractor(proteinResolver, Lists.newArrayList(), HmfGenePanelSupplier.allGenesMap37());
        return viccExtractor.extractFromViccEntries(viccEntries).knownHotspots();
    }

    @NotNull
    private static List<KnownHotspot> docmHotspots(@NotNull String docmTsv, @NotNull ProteinResolver proteinResolver) throws IOException {
        List<DocmEntry> entries = DocmReader.readAndCurate(docmTsv);

        return new DocmExtractor(proteinResolver).extractFromDocmEntries(entries);
    }

    @NotNull
    private static List<KnownHotspot> hartwigCohortHotspots(@NotNull String hartwigCohortTsv, @NotNull ProteinResolver proteinResolver,
            boolean addExplicitHotspots) throws IOException {
        LOGGER.info("Reading Hartwig Cohort TSV from '{}'", hartwigCohortTsv);
        List<HartwigEntry> hartwigCohortEntries = HartwigFileReader.read(hartwigCohortTsv);
        LOGGER.info(" Read {} entries", hartwigCohortEntries.size());

        HartwigExtractor hartwigExtractor = new HartwigExtractor(Knowledgebase.HARTWIG_COHORT, proteinResolver, addExplicitHotspots);

        return hartwigExtractor.extractFromHartwigEntries(hartwigCohortEntries);
    }

    @NotNull
    private static List<KnownHotspot> hartwigCuratedHotspots(@NotNull String hartwigCuratedTsv, @NotNull ProteinResolver proteinResolver,
            boolean addExplicitHotspots) throws IOException {
        LOGGER.info("Reading Hartwig Curated TSV from '{}'", hartwigCuratedTsv);
        List<HartwigEntry> hartwigCuratedEntries = HartwigFileReader.read(hartwigCuratedTsv);
        LOGGER.info(" Read {} entries", hartwigCuratedEntries.size());

        HartwigExtractor hartwigExtractor = new HartwigExtractor(Knowledgebase.HARTWIG_CURATED, proteinResolver, addExplicitHotspots);

        return hartwigExtractor.extractFromHartwigEntries(hartwigCuratedEntries);
    }
}
