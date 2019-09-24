package com.hartwig.hmftools.sage;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidence;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class SageVCF implements AutoCloseable {

    private final static String PASS = "PASS";

    private final String tumorSample;
    private final String normalSample;
    private final VariantContextWriter writer;

    public SageVCF(@NotNull final String filename, @NotNull final String normalSample, @NotNull final String tumorSample) {
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;

        final VCFHeader header = header(normalSample, tumorSample);
        writer = new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, false).build();
        writer.setHeader(header);
        writer.writeHeader(header);

    }

    public void write(@NotNull final List<VariantHotspotEvidence> evidenceList) {
        evidenceList.stream().map(this::create).forEach(writer::add);
    }

    public void write(@NotNull final VariantHotspotEvidence evidence) {
        writer.add(create(evidence));
    }

    @NotNull
    VariantContext create(@NotNull final VariantHotspotEvidence hotspotEvidence) {

        final Allele ref = Allele.create(hotspotEvidence.ref(), true);
        final Allele alt = Allele.create(hotspotEvidence.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final Genotype tumor = new GenotypeBuilder(tumorSample).AD(new int[] { hotspotEvidence.refSupport(), hotspotEvidence.altSupport() })
                .alleles(alleles)
                .make();

        final Genotype normal = new GenotypeBuilder(normalSample)//.DP(hotspotEvidence.normalReads())
                .AD(new int[] { 0, 0 }).alleles(alleles).make();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(hotspotEvidence.chromosome())
                .start(hotspotEvidence.position())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, round(hotspotEvidence.vaf()))
                .computeEndFromAlleles(alleles, (int) hotspotEvidence.position())
                .genotypes(tumor)
                .alleles(alleles);

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(hotspotEvidence.altQuality() / -10d);
        return context;
    }

    private static double round(double number) {
        double multiplier = Math.pow(10, 3);
        return Math.round(number * multiplier) / multiplier;
    }

    @NotNull
    private static VCFHeader header(@NotNull final String normalSample, @NotNull final String tumorSample) {
        VCFHeader header = new VCFHeader(Collections.emptySet(), Lists.newArrayList(normalSample, tumorSample));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine((VCFConstants.ALLELE_FREQUENCY_KEY)));

        return header;
    }

    @Override
    public void close() {
        writer.close();
    }

}
