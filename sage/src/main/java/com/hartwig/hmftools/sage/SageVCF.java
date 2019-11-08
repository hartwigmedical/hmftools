package com.hartwig.hmftools.sage;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.context.AltContext;
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

    public final static String PASS = "PASS";
    public final static String MERGE_FILTER = "merge";

    private final static String READ_CONTEXT = "RC";
    private final static String READ_CONTEXT_DESCRIPTION = "Read context";
    private final static String READ_CONTEXT_COUNT = "RC_CNT";
    private final static String READ_CONTEXT_COUNT_DESCRIPTION =
            "Read context counts [Full, Partial, Realigned, Shortened, Lengthened, Coverage]";
    private static final String READ_CONTEXT_REPEAT_COUNT = "RC_REPC";
    private static final String READ_CONTEXT_REPEAT_COUNT_DESCRIPTION = "Repeat count at read context";
    private static final String READ_CONTEXT_REPEAT_SEQUENCE = "RC_REPS";
    private static final String READ_CONTEXT_REPEAT_SEQUENCE_DESCRIPTION = "Repeat sequence at read context";
    private static final String READ_CONTEXT_MICRO_HOMOLOGY = "RC_MH";
    private static final String READ_CONTEXT_MICRO_HOMOLOGY_DESCRIPTION = "Micro-homology at read context";
    private static final String READ_CONTEXT_QUALITY = "RC_QUAL";
    private static final String READ_CONTEXT_QUALITY_DESCRIPTION = "Read context quality [Qual, BaseQual, MapQual, JitterPenalty]";
    private static final String READ_CONTEXT_AF_DESCRIPTION =
            "Allelic frequency calculated from read context counts as (Full + Partial + Realigned) / Coverage";

    private static final String READ_CONTEXT_DISTANCE = "RC_DIS";
    private static final String READ_CONTEXT_DISTANCE_DESCRIPTION = "Distance from read context to ref sequence";
    private static final String READ_CONTEXT_DIFFERENCE = "RC_DIF";
    private static final String READ_CONTEXT_DIFFERENCE_DESCRIPTION = "Difference between read context and ref sequence";

    private final static String TIER = "TIER";
    private final static String TIER_DESCRIPTION = "Tier: [HOTSPOT,PANEL,WIDE]";
    private final static String PHASE = "LPS";
    private final static String PHASE_DESCRIPTION = "Local Phase Set";

    private final SageConfig config;
    private final VariantContextWriter writer;
    private final SomaticRefContextEnrichment refContextEnrichment;

    SageVCF(@NotNull final IndexedFastaSequenceFile reference, @NotNull final SageConfig config) {
        this.config = config;

        writer = new VariantContextWriterBuilder().setOutputFile(config.outputFile()).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        refContextEnrichment = new SomaticRefContextEnrichment(reference, this::write);
        final VCFHeader header = refContextEnrichment.enrichHeader(header(config.reference(), config.tumor()));
        writer.writeHeader(header);
    }

    public void write(@NotNull final SageVariant entry) {
        if (shouldWrite(entry)) {
            refContextEnrichment.accept(create(entry));
        }
    }

    private void write(@NotNull final VariantContext context) {
        writer.add(context);
    }

    private boolean shouldWrite(@NotNull final SageVariant entry) {
        if (entry.isPassing()) {
            return true;
        }

        if (config.filter().hardFilter()) {
            return false;
        }

        final AltContext normal = entry.normal();
        return normal.altSupport() > config.filter().hardMaxNormalAltSupport();
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

        final ReadContextCounter firstTumorCounter = firstTumor.primaryReadContext();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(normal.chromosome())
                .start(normal.position())
                .attribute(READ_CONTEXT, firstTumorCounter.toString())
                .attribute(READ_CONTEXT_DIFFERENCE, firstTumorCounter.readContext().distanceCigar())
                .attribute(READ_CONTEXT_DISTANCE, firstTumorCounter.readContext().distance())
                .attribute(TIER, entry.tier())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, firstTumorCounter.vaf())
                .computeEndFromAlleles(alleles, (int) normal.position())
                .source("SAGE")
                .genotypes(genotypes)
                .alleles(alleles)
                .filters(entry.filters());

        if (!firstTumor.primaryReadContext().readContext().microhomology().isEmpty()) {
            builder.attribute(READ_CONTEXT_MICRO_HOMOLOGY, firstTumorCounter.readContext().microhomology());
        }

        if (firstTumor.primaryReadContext().readContext().repeatCount() > 0) {
            builder.attribute(READ_CONTEXT_REPEAT_COUNT, firstTumorCounter.readContext().repeatCount())
                    .attribute(READ_CONTEXT_REPEAT_SEQUENCE, firstTumorCounter.readContext().repeat());
        }

        if (entry.localPhaseSet() > 0) {
            builder.attribute(PHASE, entry.localPhaseSet());
        }

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(firstTumorCounter.quality() / -10d);
        if (context.isNotFiltered()) {
            context.getCommonInfo().addFilter(PASS);
        }

        return context;
    }

    @NotNull
    private Genotype createGenotype(@NotNull final List<Allele> alleles, @NotNull final AltContext evidence) {
        ReadContextCounter readContextCounter = evidence.primaryReadContext();

        return new GenotypeBuilder(evidence.sample()).DP(evidence.readDepth())
                .AD(new int[] { evidence.refSupport(), evidence.altSupport() })
                .attribute(READ_CONTEXT_QUALITY, readContextCounter.qual())
                .attribute(READ_CONTEXT_COUNT, readContextCounter.rcc())
                .alleles(alleles)
                .make();
    }

    @NotNull
    private static VCFHeader header(@NotNull final String normalSample, @NotNull final List<String> tumorSamples) {
        final List<String> allSamples = Lists.newArrayList(normalSample);
        allSamples.addAll(tumorSamples);

        VCFHeader header = new VCFHeader(Collections.emptySet(), allSamples);
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));
        header.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY,
                1,
                VCFHeaderLineType.Float,
                READ_CONTEXT_AF_DESCRIPTION));

        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_COUNT, 6, VCFHeaderLineType.Integer, READ_CONTEXT_COUNT_DESCRIPTION));
        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_QUALITY,
                4,
                VCFHeaderLineType.Integer,
                READ_CONTEXT_QUALITY_DESCRIPTION));

        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT, 1, VCFHeaderLineType.String, READ_CONTEXT_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_DIFFERENCE,
                1,
                VCFHeaderLineType.String,
                READ_CONTEXT_DIFFERENCE_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_DISTANCE,
                1,
                VCFHeaderLineType.Integer,
                READ_CONTEXT_DISTANCE_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_REPEAT_COUNT,
                1,
                VCFHeaderLineType.Integer,
                READ_CONTEXT_REPEAT_COUNT_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_REPEAT_SEQUENCE,
                1,
                VCFHeaderLineType.String,
                READ_CONTEXT_REPEAT_SEQUENCE_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_MICRO_HOMOLOGY,
                1,
                VCFHeaderLineType.String,
                READ_CONTEXT_MICRO_HOMOLOGY_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(PHASE, 1, VCFHeaderLineType.Integer, PHASE_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESCRIPTION));

        header.addMetaDataLine(new VCFFilterHeaderLine(MERGE_FILTER, "Variant was merged into another variant"));
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
        writer.close();
    }

}
