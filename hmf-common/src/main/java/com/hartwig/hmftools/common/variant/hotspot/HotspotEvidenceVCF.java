package com.hartwig.hmftools.common.variant.hotspot;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.distribution.BinomialDistribution;
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

    private static final String PASS = "PASS";
    private static final String ALLELIC_FREQUENCY = "AF";
    private static final String HOTSPOT_FLAG = "HOTSPOT";
    private static final String GERMLINE_INDEL = "GERMLINE_INDEL";
    private static final String GERMLINE_HET_LIKELIHOOD = "GHBL";
    private static final String LOW_CONFIDENCE = "LOW_CONFIDENCE";

    private final VCFHeader header;
    private final String tumorSample;
    private final String normalSample;

    private final double maxNormalHetLikelihood;
    private final int minTumorReads;
    private final double minSnvVAF;
    private final double minIndelVAF;
    private final int minSnvQuality;
    private final int minIndelQuality;

    public HotspotEvidenceVCF(@NotNull final String normalSample, @NotNull final String tumorSample, final double maxNormalHetLikelihood,
            final int minTumorReads, final double minSnvVAF, final double minIndelVAF, final int minSnvQuality, final int minIndelQuality) {
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;
        this.maxNormalHetLikelihood = maxNormalHetLikelihood;
        this.minTumorReads = minTumorReads;
        this.minSnvVAF = minSnvVAF;
        this.minIndelVAF = minIndelVAF;
        this.minSnvQuality = minSnvQuality;
        this.minIndelQuality = minIndelQuality;

        this.header = header(normalSample, tumorSample);
    }

    public void write(@NotNull final String filename, @NotNull final List<HotspotEvidence> evidenceList) {
        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, false).build();
        writer.setHeader(header);
        writer.writeHeader(header);

        final ListMultimap<GenomePosition, HotspotEvidence> evidenceMap = Multimaps.index(evidenceList, GenomePositions::create);
        for (GenomePosition site : evidenceMap.keySet()) {
            final List<HotspotEvidence> evidence = evidenceMap.get(site);
            final VariantContext context = create(evidence);
            writer.add(context);
        }

        writer.close();
    }

    private boolean lowConfidence(@NotNull HotspotEvidence hotspotEvidence) {
        if (hotspotEvidence.normalAltCount() > 1 || hotspotEvidence.tumorAltCount() < minTumorReads) {
            return true;
        }

        if (hotspotEvidence.isIndel()) {
            return hotspotEvidence.qualityScore() < minIndelQuality || Doubles.lessThan(hotspotEvidence.vaf(), minIndelVAF);
        }

        return hotspotEvidence.qualityScore() < minSnvQuality || Doubles.lessThan(hotspotEvidence.vaf(), minSnvVAF);
    }

    private static boolean germlineIndel(@NotNull HotspotEvidence hotspotEvidence) {
        return hotspotEvidence.isIndel() && hotspotEvidence.normalIndelCount() > 0;
    }

    @VisibleForTesting
    @NotNull
    VariantContext create(@NotNull final Collection<HotspotEvidence> evidence) {
        assert (!evidence.isEmpty());

        final List<VariantContext> contexts =
                evidence.stream().map(this::create).sorted(HotspotEvidenceVCF::compareEvidence).collect(Collectors.toList());

        return contexts.get(0);
    }

    @NotNull
    private VariantContext create(@NotNull final HotspotEvidence hotspotEvidence) {
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
                .attribute(ALLELIC_FREQUENCY, round(hotspotEvidence.vaf()))
                .computeEndFromAlleles(alleles, (int) hotspotEvidence.position())
                .source(hotspotEvidence.type().toString())
                .genotypes(tumor, normal)
                .alleles(alleles);

        boolean lowConfidence = lowConfidence(hotspotEvidence);
        if (hotspotEvidence.normalAltCount() == 1) {
            double normalHetLikelihood = heterozygousLikelihood(hotspotEvidence.normalReads());
            builder.attribute(GERMLINE_HET_LIKELIHOOD, round(normalHetLikelihood));
            lowConfidence |= Doubles.greaterThan(normalHetLikelihood, maxNormalHetLikelihood);
        }

        if (lowConfidence) {
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

    private static boolean isNotFiltered(@NotNull final VariantContext context) {
        return context.getFilters().contains(PASS);
    }

    private static int compareEvidence(@NotNull final VariantContext o1, @NotNull final VariantContext o2) {
        int passCompare = -Boolean.compare(isNotFiltered(o1), isNotFiltered(o2));
        return passCompare == 0 ? -Double.compare(o1.getPhredScaledQual(), o2.getPhredScaledQual()) : passCompare;
    }

    static double heterozygousLikelihood(int readDepth) {
        return new BinomialDistribution(readDepth, 0.5).cumulativeProbability(1);
    }

    private static double round(double number) {
        double multiplier = Math.pow(10, 3);
        return Math.round(number * multiplier) / multiplier;
    }

    @NotNull
    private static VCFHeader header(@NotNull final String normalSample, @NotNull final String tumorSample) {
        VCFHeader header = new VCFHeader(Collections.emptySet(), Lists.newArrayList(normalSample, tumorSample));
        header.addMetaDataLine(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read Depth"));
        header.addMetaDataLine(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Allelic Depth"));

        header.addMetaDataLine(new VCFInfoHeaderLine(ALLELIC_FREQUENCY,
                VCFHeaderLineCount.A,
                VCFHeaderLineType.Float,
                "Allelic Frequency"));
        header.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT_FLAG, 1, VCFHeaderLineType.String, "Hotspot Type: known, inframe"));
        header.addMetaDataLine(new VCFInfoHeaderLine(GERMLINE_HET_LIKELIHOOD,
                1,
                VCFHeaderLineType.Float,
                "Germline Heterozygous Binomial Likelihood. Only applies when germline alt support = 1"));

        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));
        header.addMetaDataLine(new VCFFilterHeaderLine(LOW_CONFIDENCE,
                "Set if excessive germline reads or insufficient quality or tumor reads"));
        header.addMetaDataLine(new VCFFilterHeaderLine(GERMLINE_INDEL, "Set if indel has any germline indels at that site"));

        return header;
    }
}
