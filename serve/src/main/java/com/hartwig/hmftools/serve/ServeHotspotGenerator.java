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
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.hotspot.HotspotGenerator;
import com.hartwig.hmftools.serve.util.ProteinKeyFormatter;
import com.hartwig.hmftools.serve.vicc.ViccExtractionResult;
import com.hartwig.hmftools.serve.vicc.ViccExtractor;
import com.hartwig.hmftools.serve.vicc.ViccExtractorFactory;
import com.hartwig.hmftools.serve.vicc.ViccReader;
import com.hartwig.hmftools.serve.vicc.ViccUtil;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
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

public class ServeHotspotGenerator {

    private static final Logger LOGGER = LogManager.getLogger(ServeHotspotGenerator.class);

    private static final Integer MAX_VICC_ENTRIES = null;

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String viccJsonPath;
        String refGenomeFastaFile;
        boolean generateHotspots;
        String hotspotVcf = null;

        if (hostname.toLowerCase().contains("datastore")) {
            viccJsonPath = "/data/common/dbs/serve/vicc/all.json";
            refGenomeFastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            generateHotspots = true;
            hotspotVcf = System.getProperty("user.home") + "/tmp/hotspotsServe.vcf";
        } else {
            viccJsonPath = System.getProperty("user.home") + "/hmf/projects/serve/vicc/all.json";
            refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";

            generateHotspots = false;
        }

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        LOGGER.debug("Configured '{}' as the VICC json path", viccJsonPath);
        LOGGER.debug("Configured '{}' as the reference fasta path", refGenomeFastaFile);
        LOGGER.debug("Configured '{}' as the hotspot output VCF", hotspotVcf);
        LOGGER.debug("Configured '{}' for generating hotspots yes/no", generateHotspots);

        List<ViccSource> sources = Lists.newArrayList(ViccSource.CIVIC, ViccSource.JAX, ViccSource.ONCOKB, ViccSource.CGI);
        ViccQuerySelection querySelection =
                ImmutableViccQuerySelection.builder().sourcesToFilterOn(sources).maxEntriesToInclude(MAX_VICC_ENTRIES).build();
        List<ViccEntry> viccEntries = ViccReader.readAndCurate(viccJsonPath, querySelection);

        HotspotGenerator hotspotGenerator =
                generateHotspots ? HotspotGenerator.transvarWithRefGenome(refGenomeVersion, refGenomeFastaFile) : HotspotGenerator.dummy();

        ViccExtractor viccExtractor = ViccExtractorFactory.buildViccExtractor(hotspotGenerator);
        Map<ViccEntry, ViccExtractionResult> resultsPerEntry = viccExtractor.extractFromViccEntries(viccEntries);
        Map<VariantHotspot, HotspotAnnotation> hotspotMap = ViccUtil.convertToHotspotMap(resultsPerEntry);

        if (generateHotspots && hotspotVcf != null) {
            writeHotspots(hotspotVcf, hotspotMap);

            Set<String> unresolvedProteinAnnotations = hotspotGenerator.unresolvedProteinAnnotations();
            if (!unresolvedProteinAnnotations.isEmpty()) {
                LOGGER.info("Hotspot generator could not resolve {} protein annotations", unresolvedProteinAnnotations.size());
                for (String unresolvedProteinAnnotation : unresolvedProteinAnnotations) {
                    LOGGER.warn("Hotspot generator could not resolve protein annotation '{}'", unresolvedProteinAnnotation);
                }
            } else {
                LOGGER.info("Hotspot generator could resolve hotspots for every protein annotation");
            }
        }
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
}
