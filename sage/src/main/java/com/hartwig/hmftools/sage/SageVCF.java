package com.hartwig.hmftools.sage;

import static htsjdk.variant.vcf.VCFHeaderLineCount.UNBOUNDED;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.sage.evidence.SampleEvidence;
import com.hartwig.hmftools.sage.evidence.VariantEvidence;

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
    private final static String READ_CONTEXT_COUNT_OTHER = "RCCO";

    private final VariantContextWriter writer;
    private final SomaticRefContextEnrichment refContextEnrichment;

    public SageVCF(@NotNull final IndexedFastaSequenceFile reference, @NotNull final String filename, @NotNull final String normalSample,
            @NotNull final List<String> tumorSamples) {

        writer = new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, false).build();
        refContextEnrichment = new SomaticRefContextEnrichment(reference, writer::add);

        final VCFHeader header = refContextEnrichment.enrichHeader(header(normalSample, tumorSamples));
        writer.writeHeader(header);
    }

    public void write(@NotNull final VariantEvidence evidence) {
        refContextEnrichment.accept(create(evidence));
    }

    @NotNull
    private Genotype createGenotype(@NotNull final List<Allele> alleles, @NotNull final SampleEvidence evidence) {
        return new GenotypeBuilder(evidence.sample()).DP(evidence.readDepth())
                .AD(new int[] { evidence.refSupport(), evidence.altSupport() })
                .attribute("SDP", evidence.subprimeReadDepth())
                .alleles(alleles)
                .make();
    }

    @NotNull
    private VariantContext create(@NotNull final VariantEvidence evidence) {

        final Allele ref = Allele.create(evidence.ref(), true);
        final Allele alt = Allele.create(evidence.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final List<Genotype> genotypes = Lists.newArrayList(createGenotype(alleles, evidence.normalEvidence()));
        evidence.tumorEvidence().stream().map(x -> createGenotype(alleles, x)).forEach(genotypes::add);

        final VariantContextBuilder builder = new VariantContextBuilder().chr(evidence.chromosome()).start(evidence.position())
                //                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, round(tumorEvidence.vaf()))
                //                .attribute("MAP_Q", tumorEvidence.altMapQuality())
                //                .attribute("MAP_BASE_Q", tumorEvidence.altMinQuality())
                //                .attribute("AVG_DISTANCE_RECORD", tumorEvidence.avgAltDistanceFromRecordStart())
                //                .attribute("AVG_DISTANCE_ALIGNMENT", tumorEvidence.avgAltMinDistanceFromAlignment())
                //                .attribute(READ_CONTEXT, tumorEvidence.readContext())
                //                .attribute(READ_CONTEXT_COUNT, tumorEvidence.readContextCount())
                //                .attribute(READ_CONTEXT_COUNT_OTHER, tumorEvidence.readContextCountOther())
                                .computeEndFromAlleles(alleles, (int) evidence.position())
                .source("SAGE").genotypes(genotypes).alleles(alleles);

        final VariantContext context = builder.make();
        //        context.getCommonInfo().setLog10PError(evidence.altQuality() / -10d);
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
        header.addMetaDataLine(new VCFInfoHeaderLine("MAP_Q", UNBOUNDED, VCFHeaderLineType.Float, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("MAP_BASE_Q", UNBOUNDED, VCFHeaderLineType.Float, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("AVG_DISTANCE_RECORD", UNBOUNDED, VCFHeaderLineType.Float, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("AVG_DISTANCE_ALIGNMENT", UNBOUNDED, VCFHeaderLineType.Float, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT, 1, VCFHeaderLineType.String, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_COUNT, 1, VCFHeaderLineType.Integer, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_COUNT_OTHER, 1, VCFHeaderLineType.Integer, "TODO"));

        return header;
    }

    @Override
    public void close() {
        writer.close();
    }

}
