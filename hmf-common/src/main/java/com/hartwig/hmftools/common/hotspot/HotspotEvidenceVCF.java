package com.hartwig.hmftools.common.hotspot;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;

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

    private final static int MIN_HOTSPOT_QUALITY = 100;
    private final static int MIN_INFRAME_QUALITY = 150;
    private final static int MIN_TUMOR_EVIDENCE = 2;
    private final static String LOW_CONFIDENCE = "LOW_CONFIDENCE";

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
        header.addMetaDataLine(new VCFInfoHeaderLine("HT",
                VCFHeaderLineCount.A,
                VCFHeaderLineType.String,
                "Hotspot Type: INFRAME, SNV, MNV, INSERT or DELETE"));
        header.addMetaDataLine(new VCFFilterHeaderLine(LOW_CONFIDENCE,
                "Set if not true: AD[NormalAlt] = 0 && AD[TumorAlt] >= " + MIN_TUMOR_EVIDENCE + " && QUAL[Hotspot|Inframe] >= "
                        + MIN_HOTSPOT_QUALITY + "|" + MIN_INFRAME_QUALITY));
    }

    public void write(@NotNull final String filename, @NotNull final List<HotspotEvidence> evidence) {
        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, false).build();
        writer.setHeader(header);
        writer.writeHeader(header);

        final ListMultimap<RefPosition, HotspotEvidence> evidenceMultimap = Multimaps.index(evidence, HotspotEvidenceVCF::create);
        for (RefPosition ref : evidenceMultimap.keySet()) {
            final List<HotspotEvidence> refEvidence = evidenceMultimap.get(ref);
            final VariantContext context = create(refEvidence);
            writer.add(context);
        }

        writer.close();
    }

    private static boolean lowConfidence(@NotNull HotspotEvidence hotspotEvidence) {
        return hotspotEvidence.type() == HotspotEvidenceType.INFRAME && hotspotEvidence.qualityScore() < MIN_INFRAME_QUALITY
                || hotspotEvidence.type() != HotspotEvidenceType.INFRAME && hotspotEvidence.qualityScore() < MIN_HOTSPOT_QUALITY
                || hotspotEvidence.normalAltCount() > 0 || hotspotEvidence.tumorAltCount() < MIN_TUMOR_EVIDENCE;
    }

    @VisibleForTesting
    @NotNull
    VariantContext create(@NotNull final Collection<HotspotEvidence> evidence) {
        assert (!evidence.isEmpty());

        List<HotspotEvidence> sortedEvidence = Lists.newArrayList(evidence);

        sortedEvidence.sort(HotspotEvidenceVCF::compareEvidence);

        final HotspotEvidence primaryEvidence = sortedEvidence.get(0);
        final List<Allele> alleles = Lists.newArrayList(Allele.create(primaryEvidence.ref(), true));
        sortedEvidence.stream().map(x -> Allele.create(x.alt(), false)).forEach(alleles::add);

        final List<Integer> tumorAD = Lists.newArrayList(primaryEvidence.tumorRefCount());
        sortedEvidence.stream().mapToInt(HotspotEvidence::tumorAltCount).forEach(tumorAD::add);

        final List<Integer> normalAD = Lists.newArrayList(primaryEvidence.normalRefCount());
        sortedEvidence.stream().mapToInt(HotspotEvidence::normalAltCount).forEach(normalAD::add);

        final Genotype tumor = new GenotypeBuilder(tumorSample).DP(primaryEvidence.tumorReads())
                .AD(tumorAD.stream().mapToInt(x -> x).toArray())
                .alleles(alleles)
                .make();

        final Genotype normal = new GenotypeBuilder(normalSample).DP(primaryEvidence.normalReads())
                .AD(normalAD.stream().mapToInt(x -> x).toArray())
                .alleles(alleles)
                .make();

        final StringJoiner htJoiner = new StringJoiner(",");
        sortedEvidence.stream().map(x -> x.type().toString()).forEach(htJoiner::add);

        final VariantContextBuilder builder = new VariantContextBuilder().chr(primaryEvidence.chromosome())
                .start(primaryEvidence.position())
                .computeEndFromAlleles(alleles, (int) primaryEvidence.position())
                .attribute("HT", htJoiner.toString())
                .source(primaryEvidence.type().toString())
                .genotypes(tumor, normal)
                .alleles(alleles);

        if (lowConfidence(primaryEvidence)) {
            builder.filter(LOW_CONFIDENCE);
        }

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(primaryEvidence.qualityScore() / -10d);
        return context;
    }

    @NotNull
    private static RefPosition create(@NotNull final HotspotEvidence evidence) {
        return ImmutableRefPosition.builder().from(evidence).ref(evidence.ref()).build();
    }

    private static int compareEvidence(@NotNull final HotspotEvidence o1, @NotNull final HotspotEvidence o2) {
        int normalEvidence = Integer.compare(o1.normalAltCount(), o2.normalAltCount());
        return normalEvidence == 0 ? -Integer.compare(o1.tumorAltCount(), o2.tumorAltCount()) : normalEvidence;
    }

}
