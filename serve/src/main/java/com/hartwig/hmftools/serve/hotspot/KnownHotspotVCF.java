package com.hartwig.hmftools.serve.hotspot;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class KnownHotspotVCF {

    private static final Logger LOGGER = LogManager.getLogger(KnownHotspotVCF.class);

    private static final String INFO_SOURCES = "sources";
    private static final String INFO_INPUT = "input";

    private KnownHotspotVCF() {
    }

    public static void write(@NotNull String hotspotVcf, @NotNull List<KnownHotspot> hotspots) {
        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(hotspotVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        header.addMetaDataLine(new VCFInfoHeaderLine(INFO_INPUT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "input"));
        header.addMetaDataLine(new VCFInfoHeaderLine(INFO_SOURCES,
                VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.String,
                "sources [" + uniqueSourcesString(hotspots) + "]"));

        writer.writeHeader(header);

        for (KnownHotspot hotspot : sort(hotspots)) {
            List<Allele> hotspotAlleles = buildAlleles(hotspot);

            VariantContext variantContext = new VariantContextBuilder().noGenotypes()
                    .source("SERVE")
                    .chr(hotspot.chromosome())
                    .start(hotspot.position())
                    .alleles(hotspotAlleles)
                    .computeEndFromAlleles(hotspotAlleles, (int) hotspot.position())
                    .attribute(INFO_SOURCES, Knowledgebase.commaSeparatedSourceString(hotspot.sources()))
                    .attribute(INFO_INPUT,
                            ProteinKeyFormatter.toProteinKey(hotspot.gene(), hotspot.transcript(), hotspot.proteinAnnotation()))
                    .make();

            LOGGER.debug(" Writing variant '{}'", variantContext);
            writer.add(variantContext);

        }
        writer.close();
    }

    @NotNull
    private static List<KnownHotspot> sort(@NotNull List<KnownHotspot> hotspots) {
        // Need to make a copy since the input list may be immutable and cannot be sorted!
        List<KnownHotspot> sorted = Lists.newArrayList(hotspots);
        sorted.sort(new VariantHotspotComparator());

        return sorted;
    }

    @VisibleForTesting
    @NotNull
    static String uniqueSourcesString(@NotNull List<KnownHotspot> hotspots) {
        Set<Knowledgebase> sources = Sets.newHashSet();
        for (KnownHotspot hotspot : hotspots) {
            sources.addAll(hotspot.sources());
        }
        return Knowledgebase.commaSeparatedSourceString(sources);
    }

    @NotNull
    private static List<Allele> buildAlleles(@NotNull VariantHotspot hotspot) {
        Allele ref = Allele.create(hotspot.ref(), true);
        Allele alt = Allele.create(hotspot.alt(), false);

        return Lists.newArrayList(ref, alt);
    }
}
