package com.hartwig.hmftools.serve.extraction.hotspot;

import java.io.File;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.serve.extraction.util.KeyFormatter;
import com.hartwig.hmftools.serve.extraction.util.VCFWriterFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

public final class KnownHotspotFile {

    private static final Logger LOGGER = LogManager.getLogger(KnownHotspotFile.class);
    private static final String KNOWN_HOTSPOT_VCF = "KnownHotspots.SERVE.vcf.gz";

    private KnownHotspotFile() {
    }

    @NotNull
    public static String knownHotspotVcfPath(@NotNull String outputDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(outputDir + File.separator + KNOWN_HOTSPOT_VCF);
    }

    public static void write(@NotNull String hotspotVcf, @NotNull IndexedFastaSequenceFile refSequence,
            @NotNull Iterable<KnownHotspot> hotspots) {
        VariantContextWriter writer =
                VCFWriterFactory.openIndexedVCFWriter(hotspotVcf, refSequence, uniqueSourcesString(hotspots));

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
