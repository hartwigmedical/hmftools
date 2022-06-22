package com.hartwig.hmftools.serve.extraction.hotspot;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.serve.extraction.util.KeyFormatter;
import com.hartwig.hmftools.serve.extraction.util.VCFWriterFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public final class KnownHotspotFile {

    private static final Logger LOGGER = LogManager.getLogger(KnownHotspotFile.class);
    private static final String KNOWN_HOTSPOT_VCF = "KnownHotspots.SERVE.vcf.gz";

    private KnownHotspotFile() {
    }

    @NotNull
    public static String knownHotspotVcfPath(@NotNull String outputDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(outputDir + File.separator + KNOWN_HOTSPOT_VCF);
    }

    public static List<KnownHotspot> read(@NotNull String vcfFile) throws IOException {
        final List<KnownHotspot> result = Lists.newArrayList();
        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFile,
                new VCFCodec(),
                false);) {
            for (VariantContext hotspot : reader.iterator()) {
                VariantImpact impact = VariantImpactSerialiser.fromVariantContext(hotspot);
                StringJoiner joiner = new StringJoiner("|");
                Object input = hotspot.getAttribute(VCFWriterFactory.INPUT_FIELD);
                if (input != null) {
                    joiner.add(input.toString());
                }

                List<String> sources = hotspot.getAttributeAsStringList(VCFWriterFactory.SOURCES_FIELD, Strings.EMPTY);
                if (!sources.isEmpty()) {
                    StringJoiner sourceJoiner = new StringJoiner(",");
                    for (String source : sources) {
                        sourceJoiner.add(source);
                    }
                    joiner.add(sourceJoiner.toString());
                } else {
                    LOGGER.warn("No sources found on {}", hotspot);
                }

                //TODO: add info
                result.add(ImmutableKnownHotspot.builder()
                        .chromosome(hotspot.getContig())
                        .gene(impact.CanonicalGeneName)
                        .proteinAnnotation(impact.CanonicalHgvsProtein)
                        .position(hotspot.getStart())
                        .ref(hotspot.getAlleles().get(0).getBaseString())
                        .alt(hotspot.getAlleles().get(0).getBaseString())
                        .build());
            }
        }

        return result;

    }

    public static void write(@NotNull String hotspotVcf, @NotNull IndexedFastaSequenceFile refSequence,
            @NotNull Iterable<KnownHotspot> hotspots) {
        VariantContextWriter writer = VCFWriterFactory.openIndexedVCFWriter(hotspotVcf, refSequence, uniqueSourcesString(hotspots));

        for (KnownHotspot hotspot : sort(hotspots)) {
            List<Allele> hotspotAlleles = buildAlleles(hotspot);

            VariantContext variant = new VariantContextBuilder().noGenotypes()
                    .source("SERVE")
                    .chr(hotspot.chromosome())
                    .start(hotspot.position())
                    .alleles(hotspotAlleles)
                    .computeEndFromAlleles(hotspotAlleles, hotspot.position())
                    .attribute(VCFWriterFactory.INPUT_FIELD,
                            KeyFormatter.toProteinKey(hotspot.gene(), hotspot.transcript(), hotspot.proteinAnnotation()))
                    .attribute(VCFWriterFactory.SOURCES_FIELD, Knowledgebase.toCommaSeparatedSourceString(hotspot.sources()))
                    .make();

            LOGGER.debug(" Writing variant '{}'", variant);
            writer.add(variant);
        }

        writer.close();
    }

    @NotNull
    private static List<KnownHotspot> sort(@NotNull Iterable<KnownHotspot> hotspots) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<KnownHotspot> sorted = Lists.newArrayList(hotspots);
        sorted.sort(new VariantHotspotComparator());

        return sorted;
    }

    @VisibleForTesting
    @NotNull
    static String uniqueSourcesString(@NotNull Iterable<KnownHotspot> hotspots) {
        Set<Knowledgebase> sources = Sets.newHashSet();
        for (KnownHotspot hotspot : hotspots) {
            sources.addAll(hotspot.sources());
        }
        return Knowledgebase.toCommaSeparatedSourceString(sources);
    }

    @NotNull
    private static List<Allele> buildAlleles(@NotNull VariantHotspot hotspot) {
        Allele ref = Allele.create(hotspot.ref(), true);
        Allele alt = Allele.create(hotspot.alt(), false);

        return Lists.newArrayList(ref, alt);
    }
}
