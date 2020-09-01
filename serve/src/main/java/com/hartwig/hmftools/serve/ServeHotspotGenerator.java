package com.hartwig.hmftools.serve;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.docm.DocmEntry;
import com.hartwig.hmftools.serve.docm.DocmExtractor;
import com.hartwig.hmftools.serve.docm.DocmFileReader;
import com.hartwig.hmftools.serve.docm.curation.DocmCurator;
import com.hartwig.hmftools.serve.hartwig.HartwigEntry;
import com.hartwig.hmftools.serve.hartwig.HartwigExtractor;
import com.hartwig.hmftools.serve.hartwig.cohort.HartwigCohortEntry;
import com.hartwig.hmftools.serve.hartwig.cohort.HartwigCohortFileReader;
import com.hartwig.hmftools.serve.hartwig.curated.HartwigCuratedEntry;
import com.hartwig.hmftools.serve.hartwig.curated.HartwigCuratedFileReader;
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.hotspot.ProteinKeyFormatter;
import com.hartwig.hmftools.serve.hotspot.ProteinToHotspotConverter;
import com.hartwig.hmftools.serve.vicc.ViccExtractionResult;
import com.hartwig.hmftools.serve.vicc.ViccExtractor;
import com.hartwig.hmftools.serve.vicc.ViccExtractorFactory;
import com.hartwig.hmftools.serve.vicc.ViccReader;
import com.hartwig.hmftools.serve.vicc.ViccUtil;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

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

public class ServeHotspotGenerator {

    private static final Logger LOGGER = LogManager.getLogger(ServeHotspotGenerator.class);

    private static final Set<ViccSource> VICC_SOURCES_TO_INCLUDE =
            Sets.newHashSet(ViccSource.CIVIC, ViccSource.JAX, ViccSource.ONCOKB, ViccSource.CGI);
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
            serveSourceDir = "/data/common/dbs/serve";

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

        ProteinToHotspotConverter proteinToHotspotConverter = generateHotspots
                ? ProteinToHotspotConverter.transvarWithRefGenome(refGenomeVersion, refGenomeFastaFile)
                : ProteinToHotspotConverter.dummy();

        Map<VariantHotspot, HotspotAnnotation> viccHotspotMap = viccHotspotMap(viccJson, proteinToHotspotConverter);
        Map<VariantHotspot, HotspotAnnotation> docmHotspotMap = docmHotspotMap(docmTsv, proteinToHotspotConverter);
        Map<VariantHotspot, HotspotAnnotation> hartwigCohortMap =
                hartwigCohortMap(hartwigCohortTsv, proteinToHotspotConverter, generateHotspots);
        Map<VariantHotspot, HotspotAnnotation> hartwigCuratedMap =
                hartwigCuratedMap(hartwigCuratedTsv, proteinToHotspotConverter, generateHotspots);

        LOGGER.info("Merging {} VICC hotspots with {} DoCM hotspots and {} Hartwig Cohort hotspots and {} Hartwig Curated hotspots",
                viccHotspotMap.size(),
                docmHotspotMap.size(),
                hartwigCohortMap.size(),
                hartwigCuratedMap.size());

        Map<VariantHotspot, HotspotAnnotation> mergedMap =
                HotspotFunctions.mergeHotspots(Lists.newArrayList(viccHotspotMap, docmHotspotMap, hartwigCohortMap, hartwigCuratedMap));

        if (generateHotspots && hotspotVcf != null) {
            writeHotspots(hotspotVcf, mergedMap);

            Set<String> unresolvedProteinAnnotations = proteinToHotspotConverter.unresolvedProteinAnnotations();
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
    private static Map<VariantHotspot, HotspotAnnotation> viccHotspotMap(@NotNull String viccJson,
            @NotNull ProteinToHotspotConverter proteinToHotspotConverter) throws IOException {
        List<ViccEntry> viccEntries = ViccReader.readAndCurate(viccJson, VICC_SOURCES_TO_INCLUDE, MAX_VICC_ENTRIES);
        ViccExtractor viccExtractor = ViccExtractorFactory.buildViccExtractor(proteinToHotspotConverter);
        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = viccExtractor.extractFromViccEntries(viccEntries);

        return ViccUtil.convertToHotspotMap(resultsPerEntry);
    }

    @NotNull
    private static Map<VariantHotspot, HotspotAnnotation> docmHotspotMap(@NotNull String docmTsv,
            @NotNull ProteinToHotspotConverter proteinToHotspotConverter) throws IOException {
        LOGGER.info("Reading DoCM TSV from '{}'", docmTsv);
        List<DocmEntry> docmEntries = DocmFileReader.readDcomFile(docmTsv);
        LOGGER.info(" Read {} entries", docmEntries.size());

        DocmCurator curator = new DocmCurator();
        LOGGER.info("Curating {} DoCM entries", docmEntries.size());
        List<DocmEntry> curatedDocmEntries = curator.curate(docmEntries);
        LOGGER.info(" Finished DoCM curation. {} entries remaining after curation. {} entries have been removed",
                curatedDocmEntries.size(),
                docmEntries.size() - curatedDocmEntries.size());
        curator.reportUnusedBlacklistEntries();

        DocmExtractor docmExtractor = new DocmExtractor(proteinToHotspotConverter);
        Map<DocmEntry, List<VariantHotspot>> docmHotspotsPerEntry = docmExtractor.extractFromDocmEntries(curatedDocmEntries);

        return HotspotFunctions.convertHotspots("docm", docmHotspotsPerEntry);
    }

    @NotNull
    private static Map<VariantHotspot, HotspotAnnotation> hartwigCohortMap(@NotNull String hartwigCohortTsv,
            @NotNull ProteinToHotspotConverter proteinToHotspotConverter, boolean addImpliedHotspots) throws IOException {
        LOGGER.info("Reading Hartwig Cohort TSV from '{}'", hartwigCohortTsv);
        List<HartwigCohortEntry> hartwigCohortEntries = HartwigCohortFileReader.readCohortFile(hartwigCohortTsv);
        LOGGER.info(" Read {} entries", hartwigCohortEntries.size());

        HartwigExtractor hartwigExtractor = new HartwigExtractor(proteinToHotspotConverter, addImpliedHotspots);
        Map<HartwigEntry, List<VariantHotspot>> cohortHotspotsPerEntry = hartwigExtractor.extractFromHartwigEntries(hartwigCohortEntries);

        return HotspotFunctions.convertHotspots("hartwig_cohort", cohortHotspotsPerEntry);
    }

    @NotNull
    private static Map<VariantHotspot, HotspotAnnotation> hartwigCuratedMap(@NotNull String hartwigCuratedTsv,
            @NotNull ProteinToHotspotConverter proteinToHotspotConverter, boolean addImpliedHotspots) throws IOException {
        LOGGER.info("Reading Hartwig Curated TSV from '{}'", hartwigCuratedTsv);
        List<HartwigCuratedEntry> hartwigCuratedEntries = HartwigCuratedFileReader.readCuratedFile(hartwigCuratedTsv);
        LOGGER.info(" Read {} entries", hartwigCuratedEntries.size());

        HartwigExtractor hartwigExtractor = new HartwigExtractor(proteinToHotspotConverter, addImpliedHotspots);
        Map<HartwigEntry, List<VariantHotspot>> curatedHotspotsPerEntry = hartwigExtractor.extractFromHartwigEntries(hartwigCuratedEntries);

        return HotspotFunctions.convertHotspots("hartwig_curated", curatedHotspotsPerEntry);
    }

    private static void writeHotspots(@NotNull String hotspotVcf, @NotNull Map<VariantHotspot, HotspotAnnotation> hotspotMap) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(hotspotVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        writer.writeHeader(header);

        for (Map.Entry<VariantHotspot, HotspotAnnotation> entry : hotspotMap.entrySet()) {
            VariantHotspot hotspot = entry.getKey();
            HotspotAnnotation annotation = entry.getValue();
            List<Allele> hotspotAlleles = buildAlleles(hotspot);

            VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                    .source("SERVE")
                    .chr(hotspot.chromosome())
                    .start(hotspot.position())
                    .alleles(hotspotAlleles)
                    .computeEndFromAlleles(hotspotAlleles, (int) hotspot.position())
                    .attribute("sources", buildSourcesString(annotation.sources()))
                    .attribute("input", ProteinKeyFormatter.toProteinKey(annotation))
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
}
