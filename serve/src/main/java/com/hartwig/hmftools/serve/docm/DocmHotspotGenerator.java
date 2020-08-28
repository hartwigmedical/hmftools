package com.hartwig.hmftools.serve.docm;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.hotspot.HotspotAnnotation;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.transvar.Transvar;
import com.hartwig.hmftools.serve.util.ProteinKeyFormatter;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class DocmHotspotGenerator {

    private static final Logger LOGGER = LogManager.getLogger(DocmHotspotGenerator.class);

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String docmInputTsv;
        String refGenomeFastaFile;
        boolean generateHotspots;
        String hotspotVcf = null;

        if (hostname.toLowerCase().contains("datastore")) {
            docmInputTsv = "/data/common/dbs/docm/docm_v3.2.tsv";
            refGenomeFastaFile = "/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            generateHotspots = true;
            hotspotVcf = System.getProperty("user.home") + "/tmp/hotspotsDocm.vcf";
        } else {
            docmInputTsv = System.getProperty("user.home") + "/hmf/projects/docm/docm_v3.2.tsv";
            refGenomeFastaFile = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
            generateHotspots = false;
        }

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.HG19;
        LOGGER.debug("Configured '{}' as the DoCM input TSV path", docmInputTsv);
        LOGGER.debug("Configured '{}' as the reference fasta path", refGenomeFastaFile);
        LOGGER.debug("Configured '{}' as the hotspot output VCF", hotspotVcf);
        LOGGER.debug("Configured '{}' for generating hotspots yes/no", generateHotspots);

        ProteinResolver proteinResolver =
                generateHotspots ? Transvar.withRefGenome(refGenomeVersion, refGenomeFastaFile) : new ProteinResolver() {
                    @NotNull
                    @Override
                    public List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull final String gene,
                            @Nullable final String specificTranscript, @NotNull final String proteinAnnotation) {
                        return Lists.newArrayList();
                    }

                    @NotNull
                    @Override
                    public Set<String> unresolvedProteinAnnotations() {
                        return Sets.newHashSet();
                    }
                };

        List<DocmEntry> entries = DocmFileReader.readDcomFile(docmInputTsv);
        LOGGER.info("Read {} DoCM entries from {}", entries.size(), docmInputTsv);

        Map<DocmEntry, List<VariantHotspot>> hotspotsPerEntry = Maps.newHashMap();
        for (DocmEntry entry : entries) {
            LOGGER.debug("Generating hotspots for {}",
                    ProteinKeyFormatter.toProteinKey(entry.gene(), entry.transcript(), entry.proteinAnnotation()));
            hotspotsPerEntry.put(entry,
                    proteinResolver.extractHotspotsFromProteinAnnotation(entry.gene(), entry.transcript(), entry.proteinAnnotation()));
        }

        if (generateHotspots) {
            writeHotspots(hotspotVcf, hotspotsPerEntry);
        }
    }

    private static void writeHotspots(@NotNull String hotspotVcf, @NotNull Map<DocmEntry, List<VariantHotspot>> hotspotsPerEntry) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(hotspotVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        writer.writeHeader(header);

        for (Map.Entry<VariantHotspot, HotspotAnnotation> entry : convertAndSort(hotspotsPerEntry).entrySet()) {
            VariantHotspot hotspot = entry.getKey();
            HotspotAnnotation annotation = entry.getValue();
            List<Allele> hotspotAlleles = buildAlleles(hotspot);

            VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                    .source("DoCM")
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

    @NotNull
    private static Map<VariantHotspot, HotspotAnnotation> convertAndSort(@NotNull Map<DocmEntry, List<VariantHotspot>> hotspotsPerEntry) {
        Map<VariantHotspot, HotspotAnnotation> convertedMap = Maps.newTreeMap(new VariantHotspotComparator());
        for (Map.Entry<DocmEntry, List<VariantHotspot>> entry : hotspotsPerEntry.entrySet()) {
            DocmEntry docmEntry = entry.getKey();
            for (VariantHotspot hotspot : entry.getValue()) {
                HotspotAnnotation newAnnotation = new HotspotAnnotation(Sets.newHashSet("docm"),
                        docmEntry.gene(),
                        docmEntry.transcript(),
                        docmEntry.proteinAnnotation());

                HotspotAnnotation currentAnnotation = convertedMap.get(hotspot);
                if (currentAnnotation != null) {
                    LOGGER.warn("Annotation '{}' already found previously: '{}'", newAnnotation, currentAnnotation);
                } else {
                    convertedMap.put(hotspot, newAnnotation);
                }
            }
        }

        return convertedMap;
    }
}
