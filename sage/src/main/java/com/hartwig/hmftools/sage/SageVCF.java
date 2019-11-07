package com.hartwig.hmftools.sage;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.phase.PhasingQueue;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class SageVCF implements AutoCloseable {

    private final static String SUBPRIME_QUALITY_READ_DEPTH = "SDP";
    private final static String IMPROPER_PAIR_FLAG = "IPF";

    private final static String READ_CONTEXT = "RC";
    private final static String PASS = "PASS";
    private final static String READ_CONTEXT_COUNT = "RCC";
    private final static String TIER = "TIER";
    private final static String TIER_DESCRIPTION = "Tier: [HOTSPOT,PANEL,WIDE]";
    private final static String PHASE = "LPS";

    public static final String REPEAT_COUNT_FLAG = "RCREPC";
    public static final String REPEAT_SEQUENCE_FLAG = "RCREPS";

    private final SageConfig config;
    private final VariantContextWriter writer;
    private final SomaticRefContextEnrichment refContextEnrichment;
    private final PhasingQueue phasingQueue;

    SageVCF(@NotNull final IndexedFastaSequenceFile reference, @NotNull final SageConfig config) {
        this.config = config;

        writer = new VariantContextWriterBuilder().setOutputFile(config.outputFile()).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        refContextEnrichment = new SomaticRefContextEnrichment(reference, this::write);

        final VCFHeader header = refContextEnrichment.enrichHeader(header(config.reference(), config.tumor()));
        phasingQueue = new PhasingQueue(entry -> refContextEnrichment.accept(create(entry)));

        writer.writeHeader(header);
    }

    public void write(@NotNull final SageVariant entry) {
        final AltContext normal = entry.normal();
        if (normal.altSupport() <= config.maxNormalAltSupport()) {
            phasingQueue.accept(entry);
        }
    }

    private void write(@NotNull final VariantContext context) {
        if (!config.filter().hardFilter() || context.getFilters().contains(PASS)) {
            writer.add(context);
        }
    }

    @NotNull
    private Genotype createGenotype(@NotNull final List<Allele> alleles, @NotNull final AltContext evidence) {
        ReadContextCounter readContextCounter = evidence.primaryReadContext();

        return new GenotypeBuilder(evidence.sample()).DP(readContextCounter.coverage())
                .AD(new int[] { evidence.refSupport(), readContextCounter.support() })
                //                .attribute("SDP", evidence.subprimeReadDepth())
                .attribute("QUAL", readContextCounter.qual())
                .attribute("RCC", readContextCounter.rcc())
                .alleles(alleles)
                .make();
    }

    @NotNull
    private VariantContext create(@NotNull final SageVariant entry) {
        final AltContext normal = entry.normal();
        final List<AltContext> tumorContexts = entry.tumorAltContexts();

        assert (tumorContexts.size() >= 1);

        final AltContext firstTumor = tumorContexts.get(0);

        final Allele ref = Allele.create(normal.ref(), true);
        final Allele alt = Allele.create(normal.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);
        final Genotype normalGenotype = createGenotype(alleles, normal);

        final List<Genotype> genotypes = tumorContexts.stream().map(x -> createGenotype(alleles, x)).collect(Collectors.toList());
        genotypes.add(0, normalGenotype);

        final VariantContextBuilder builder = new VariantContextBuilder().chr(normal.chromosome())
                .start(normal.position())
                .attribute("RC", normal.primaryReadContext().toString())
                .attribute("RC_DIF", normal.primaryReadContext().readContext().distanceCigar())
                .attribute("RC_DIS", normal.primaryReadContext().readContext().distance())
                .attribute("TIER", entry.tier())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, firstTumor.primaryReadContext().vaf())
                .computeEndFromAlleles(alleles, (int) normal.position())
                .source("SAGE")
                .genotypes(genotypes)
                .alleles(alleles)
                .filters(entry.filters());

        if (!firstTumor.primaryReadContext().readContext().microhomology().isEmpty()) {
            builder.attribute("RC_MH", firstTumor.primaryReadContext().readContext().microhomology());
        }

        if (firstTumor.primaryReadContext().readContext().repeatCount() > 0) {
            builder.attribute("RC_REPC", firstTumor.primaryReadContext().readContext().repeatCount())
                    .attribute("RC_REPS", firstTumor.primaryReadContext().readContext().repeat());
        }

        if (firstTumor.phase() > 0) {
            builder.attribute(PHASE, firstTumor.phase());
        }

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(firstTumor.primaryReadContext().quality() / -10d);
        return context;
    }

    private static double round(double number) {
        double multiplier = Math.pow(10, 3);
        return Math.round(number * multiplier) / multiplier;
    }

    @NotNull
    private static VCFHeader header(@NotNull final String normalSample, @NotNull final List<String> tumorSamples) {
        final List<String> allSamples = Lists.newArrayList(normalSample);
        allSamples.addAll(tumorSamples);

        VCFHeader header = new VCFHeader(Collections.emptySet(), allSamples);
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine((VCFConstants.ALLELE_FREQUENCY_KEY)));
        //        header.addMetaDataLine(new VCFFormatHeaderLine(SUBPRIME_QUALITY_READ_DEPTH,1, VCFHeaderLineType.Integer,"Subprime quality read depth"));

        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT, 1, VCFHeaderLineType.String, "Read context"));
        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_COUNT,
                5,
                VCFHeaderLineType.Integer,
                "[Full, Partial, Realigned, Shortened, Lengthened]"));
        //        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_QUALITY,
        //                3,
        //                VCFHeaderLineType.Integer,
        //                "[ImproperPairedRead, InconsistentChromosome, ExcessInferredSize]"));

        header.addMetaDataLine(new VCFFormatHeaderLine("QUAL", 3, VCFHeaderLineType.Integer, "[BaseQual, MapQual, JitterPenalty]"));
        header.addMetaDataLine(new VCFFormatHeaderLine("DIST", 2, VCFHeaderLineType.Integer, "[AvgRecordDistance, AvgAlignmentDistance]"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RC_DIF", 1, VCFHeaderLineType.String, "Difference from ref"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RC_DIS", 1, VCFHeaderLineType.Integer, "Distance from ref"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RC_REPC", 1, VCFHeaderLineType.Integer, "Repeat count in read context"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RC_REPS", 1, VCFHeaderLineType.String, "Repeat sequence in read context"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RC_MH", 1, VCFHeaderLineType.String, "Microhomology in read context"));
        header.addMetaDataLine(new VCFInfoHeaderLine(PHASE, 1, VCFHeaderLineType.Integer, "Local phase set"));
        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESCRIPTION));

        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilterConfig.MIN_TUMOR_QUAL, "Insufficient tumor quality"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilterConfig.MIN_TUMOR_VAF, "Insufficient tumor VAF"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilterConfig.MIN_GERMLINE_DEPTH, "Insufficient germline depth"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilterConfig.MAX_GERMLINE_VAF, "Excess germline VAF"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilterConfig.MAX_GERMLINE_REL_QUAL, "Excess germline relative quality"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilterConfig.MAX_GERMLINE_REL_RCC,
                "Excess germline relative read context count"));

        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));

        return header;
    }

    @Override
    public void close() {
        phasingQueue.flush();
        writer.close();
    }

}
