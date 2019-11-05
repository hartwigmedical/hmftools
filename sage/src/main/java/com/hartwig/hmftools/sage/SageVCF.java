package com.hartwig.hmftools.sage;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.sage.config.SoftFilterConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.ReadContextCounter;

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

    private final static String PASS = "PASS";
    private final static String SUBPRIME_QUALITY_READ_DEPTH = "SDP";
    private final static String READ_CONTEXT = "RC";
    private final static String READ_CONTEXT_COUNT = "RCC";
    private final static String READ_CONTEXT_QUALITY = "RCQ";
    private final static String TIER = "TIER";
    private final static String TIER_DESCRIPTION = "Tier: [HOTSPOT,PANEL,WIDE]";

    private final SageConfig config;
    private final VariantContextWriter writer;
    private final SomaticRefContextEnrichment refContextEnrichment;
    private final PhasingQueue phasingQueue;
    private final SagePanel panel;
    private final HotspotEnrichment hotspotEnrichment;

    SageVCF(@NotNull final ListMultimap<Chromosome, VariantHotspot> hotspots,
            @NotNull final ListMultimap<Chromosome, GenomeRegion> panelRegions, @NotNull final IndexedFastaSequenceFile reference,
            @NotNull final SageConfig config) {
        this.panel = new SagePanel(panelRegions);
        this.config = config;
        this.hotspotEnrichment = new HotspotEnrichment(hotspots);

        writer = new VariantContextWriterBuilder().setOutputFile(config.outputFile()).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        refContextEnrichment = new SomaticRefContextEnrichment(reference, writer::add);

        final VCFHeader header = refContextEnrichment.enrichHeader(header(config.reference(), config.tumor()));
        phasingQueue = new PhasingQueue(entry -> refContextEnrichment.accept(create(entry)));

        writer.writeHeader(header);
    }

    public void write(@NotNull final SageEntry entry) {
        final AltContext normal = entry.normal();
        if (normal.altSupport() <= config.maxNormalAltSupport()) {
            phasingQueue.accept(entry);
        }
    }

    @NotNull
    private Genotype createGenotype(@NotNull final List<Allele> alleles, @NotNull final AltContext evidence) {
        ReadContextCounter readContextCounter = evidence.primaryReadContext();

        return new GenotypeBuilder(evidence.sample()).DP(evidence.readDepth())
                .AD(new int[] { evidence.refSupport(), evidence.altSupport() })
                .attribute("SDP", evidence.subprimeReadDepth())
                .attribute("QUAL",
                        new int[] { readContextCounter.quality(), readContextCounter.baseQuality(), readContextCounter.mapQuality() })
                .attribute("RCC", readContextCounter.rcc())
                .attribute("RCQ", readContextCounter.rcq())
                .alleles(alleles)
                .make();
    }

    @NotNull
    private VariantContext create(@NotNull final SageEntry entry) {
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

                //                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, round(tumorEvidence.vaf()))
                //                .attribute("MAP_Q", tumorEvidence.altMapQuality())
                //                .attribute("MAP_BASE_Q", tumorEvidence.altMinQuality())
                //                .attribute("AVG_DISTANCE_RECORD", tumorEvidence.avgAltDistanceFromRecordStart())
                //                .attribute("AVG_DISTANCE_ALIGNMENT", tumorEvidence.avgAltMinDistanceFromAlignment())
                //                .attribute(READ_CONTEXT, tumorEvidence.readContext())
                //                .attribute(READ_CONTEXT_COUNT, tumorEvidence.readContextCount())
                //                .attribute(READ_CONTEXT_COUNT_OTHER, tumorEvidence.readContextCountOther())
                .attribute("RC", normal.primaryReadContext().toString())
                .attribute("RDIF", normal.primaryReadContext().readContext().distanceCigar())
                .attribute("RDIS", normal.primaryReadContext().readContext().distance())
                .attribute("RCMH", firstTumor.primaryReadContext().readContext().microhomology())
                .attribute("RCR", firstTumor.primaryReadContext().readContext().repeat())
                .attribute("JITTER", firstTumor.primaryReadContext().qualityJitterPenalty())
                .computeEndFromAlleles(alleles, (int) normal.position())
                .source("SAGE")
                .genotypes(genotypes)
                .alleles(alleles);

        if (firstTumor.phase() > 0) {
            builder.attribute("PHASE", firstTumor.phase());
        }


        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(firstTumor.primaryReadContext().quality() / -10d);
        if (hotspotEnrichment.isOnHotspot(context)) {
            filter(config.filter().softHotspotFilter(), entry).forEach(x -> context.getCommonInfo().addFilter(x));
            context.getCommonInfo().putAttribute(TIER, "HOTSPOT");
        } else if (panel.inPanel(normal)) {
            filter(config.filter().softPanelFilter(), entry).forEach(x -> context.getCommonInfo().addFilter(x));
            context.getCommonInfo().putAttribute(TIER, "PANEL");
        } else {
             filter(config.filter().softWideFilter(), entry).forEach(x -> context.getCommonInfo().addFilter(x));
            context.getCommonInfo().putAttribute(TIER, "WIDE");
        }
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
        header.addMetaDataLine(new VCFFormatHeaderLine(SUBPRIME_QUALITY_READ_DEPTH,
                1,
                VCFHeaderLineType.Integer,
                "Subprime quality read depth"));

        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine((VCFConstants.ALLELE_FREQUENCY_KEY)));

        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT, 1, VCFHeaderLineType.String, "TODO"));
        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_COUNT,
                5,
                VCFHeaderLineType.Integer,
                "[Full, Partial, Realigned, J-, J+]"));
        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_QUALITY,
                3,
                VCFHeaderLineType.Integer,
                "[ImproperPairedRead, InconsistentChromosome, ExcessInferredSize]"));

        header.addMetaDataLine(new VCFFormatHeaderLine("QUAL", 3, VCFHeaderLineType.Integer, "[MinBaseMapQual, BaseQual, MapQual]"));
        header.addMetaDataLine(new VCFFormatHeaderLine("DIST", 2, VCFHeaderLineType.Integer, "[AvgRecordDistance, AvgAlignmentDistance]"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RDIF", 1, VCFHeaderLineType.String, "Difference from ref"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RDIS", 1, VCFHeaderLineType.Integer, "Distance from ref"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RCR", 1, VCFHeaderLineType.String, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("RCMH", 1, VCFHeaderLineType.String, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("PHASE", 1, VCFHeaderLineType.Integer, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("JITTER", 1, VCFHeaderLineType.Integer, "TODO"));
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

    @NotNull
    private List<String> filter(@NotNull final SoftFilterConfig config, @NotNull final SageEntry entry) {
        List<String> result = Lists.newArrayList();

        final AltContext primaryTumor = entry.primaryTumor();
        if (primaryTumor.primaryReadContext().quality() < config.minTumorQual()) {
            result.add(SoftFilterConfig.MIN_TUMOR_QUAL);
        }

        final AltContext normal = entry.normal();
        if (normal.readDepth() < config.minGermlineDepth()) {
            result.add(SoftFilterConfig.MIN_GERMLINE_DEPTH);
        }

        if (result.isEmpty()) {
            result.add(PASS);
        }

        return result;
    }

    @Override
    public void close() {
        phasingQueue.flush();
        writer.close();
    }

}
