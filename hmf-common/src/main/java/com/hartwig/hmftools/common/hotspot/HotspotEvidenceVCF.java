package com.hartwig.hmftools.common.hotspot;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class HotspotEvidenceVCF {

    private final static int MIN_TUMOR_EVIDENCE = 2;
    private final static int MIN_KNOWN_QUALITY = 100;
    private final static int MIN_INFRAME_QUALITY = 150;

    private final static String PASS = "PASS";
    private final static String HOTSPOT_FLAG = "HOTSPOT";
    private final static String GERMLINE_INDEL = "GERMLINE_INDEL";
    private final static String GERMLINE_INDEL_DESCRIPTION = "Set if inframe indel has any indels at that site.";

    private final static String LOW_CONFIDENCE = "LOW_CONFIDENCE";
    private final static String LOW_CONFIDENCE_DESCRIPTION =
            "Set if not true: AD[NormalAlt] = 0 && AD[TumorAlt] >= " + MIN_TUMOR_EVIDENCE + " && QUAL[Known|Inframe] >= "
                    + MIN_KNOWN_QUALITY + "|" + MIN_INFRAME_QUALITY;

    private final String tumorSample;
    private final String normalSample;
    private final VCFHeader header;

    public HotspotEvidenceVCF(@NotNull final String normalSample, @NotNull final String tumorSample) {
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;

        this.header = new VCFHeader(Collections.emptySet(), Lists.newArrayList(normalSample, tumorSample));
        header.addMetaDataLine(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read Depth"));
        header.addMetaDataLine(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Allelic Depth"));
        header.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT_FLAG, 1, VCFHeaderLineType.String, "Hotspot Type: known, inframe"));
        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));
        header.addMetaDataLine(new VCFFilterHeaderLine(LOW_CONFIDENCE, LOW_CONFIDENCE_DESCRIPTION));
        header.addMetaDataLine(new VCFFilterHeaderLine(GERMLINE_INDEL, GERMLINE_INDEL_DESCRIPTION));
    }

    public void write(@NotNull final String filename, @NotNull final List<HotspotEvidence> evidence) {
        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, false).build();
        writer.setHeader(header);
        writer.writeHeader(header);

        final ListMultimap<GenomePosition, HotspotEvidence> evidenceMultimap = Multimaps.index(evidence, GenomePositions::create);
        for (GenomePosition ref : evidenceMultimap.keySet()) {
            final List<HotspotEvidence> refEvidence = evidenceMultimap.get(ref);
            final VariantContext context = create(refEvidence);
            writer.add(context);
        }

        writer.close();
    }

    private static boolean lowConfidence(@NotNull HotspotEvidence hotspotEvidence) {
        return hotspotEvidence.type() == HotspotEvidenceType.INFRAME && hotspotEvidence.qualityScore() < MIN_INFRAME_QUALITY
                || hotspotEvidence.type() == HotspotEvidenceType.KNOWN && hotspotEvidence.qualityScore() < MIN_KNOWN_QUALITY
                || hotspotEvidence.normalAltCount() > 0
                || hotspotEvidence.tumorAltCount() < MIN_TUMOR_EVIDENCE;
    }

    private static boolean germlineIndel(@NotNull HotspotEvidence hotspotEvidence) {
        return hotspotEvidence.isIndel() && hotspotEvidence.normalIndelCount() > 0;
    }

    @VisibleForTesting
    @NotNull
    VariantContext create(@NotNull final Collection<HotspotEvidence> evidence) {
        assert (!evidence.isEmpty());

        List<HotspotEvidence> sortedEvidence = Lists.newArrayList(evidence);
        sortedEvidence.sort(HotspotEvidenceVCF::compareEvidence);

        final HotspotEvidence hotspotEvidence = sortedEvidence.get(0);

        final Allele ref = Allele.create(hotspotEvidence.ref(), true);
        final Allele alt = Allele.create(hotspotEvidence.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final Genotype tumor = new GenotypeBuilder(tumorSample).DP(hotspotEvidence.tumorReads())
                .AD(new int[] { hotspotEvidence.tumorRefCount(), hotspotEvidence.tumorAltCount() })
                .alleles(alleles)
                .make();

        final Genotype normal = new GenotypeBuilder(normalSample).DP(hotspotEvidence.normalReads())
                .AD(new int[] { hotspotEvidence.normalRefCount(), hotspotEvidence.normalAltCount() })
                .alleles(alleles)
                .make();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(hotspotEvidence.chromosome())
                .start(hotspotEvidence.position())
                .attribute(HOTSPOT_FLAG, hotspotEvidence.type().toString().toLowerCase())
                .computeEndFromAlleles(alleles, (int) hotspotEvidence.position())
                .source(hotspotEvidence.type().toString())
                .genotypes(tumor, normal)
                .alleles(alleles);

        if (lowConfidence(hotspotEvidence)) {
            builder.filter(LOW_CONFIDENCE);
        } else if (germlineIndel(hotspotEvidence)) {
            builder.filter(GERMLINE_INDEL);
        } else {
            builder.filter(PASS);
        }

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(hotspotEvidence.qualityScore() / -10d);
        return context;
    }

    private static int compareEvidence(@NotNull final HotspotEvidence o1, @NotNull final HotspotEvidence o2) {
        int normalEvidence = Integer.compare(o1.normalAltCount(), o2.normalAltCount());
        return normalEvidence == 0 ? -Integer.compare(o1.qualityScore(), o2.qualityScore()) : normalEvidence;
    }
}
