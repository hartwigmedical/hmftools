package com.hartwig.hmftools.sage;

import static htsjdk.variant.vcf.VCFHeaderLineCount.UNBOUNDED;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidence;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
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


    public void write(@NotNull final VariantHotspotEvidence tumorEvidence, @Nullable final VariantHotspotEvidence normalEvidence) {
        writer.add(create(tumorEvidence, normalEvidence));
    }

    @NotNull
    private VariantContext create(@NotNull final VariantHotspotEvidence tumorEvidence, @Nullable final VariantHotspotEvidence normalEvidence) {

        final Allele ref = Allele.create(tumorEvidence.ref(), true);
        final Allele alt = Allele.create(tumorEvidence.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final Genotype tumor = new GenotypeBuilder(tumorSample).DP(tumorEvidence.readDepth())
                .AD(new int[] { tumorEvidence.refSupport(), tumorEvidence.altSupport() })
                .alleles(alleles)
                .make();

        final Genotype normal;
        if (normalEvidence == null) {
            normal = new GenotypeBuilder(normalSample).DP(0).AD(new int[] { 0, 0 }).alleles(alleles).make();
        } else {
            normal = new GenotypeBuilder(normalSample).DP(normalEvidence.readDepth())
                    .AD(new int[] { normalEvidence.refSupport(), normalEvidence.altSupport() })
                    .alleles(alleles)
                    .make();
        }

        final VariantContextBuilder builder = new VariantContextBuilder().chr(tumorEvidence.chromosome())
                .start(tumorEvidence.position())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, round(tumorEvidence.vaf()))
                .attribute("MAP_Q", tumorEvidence.altMapQuality())
                .attribute("MAP_BASE_Q", tumorEvidence.altMinQuality())
                .attribute("AVG_DISTANCE_RECORD", tumorEvidence.avgAltDistanceFromRecordStart())
                .attribute("AVG_DISTANCE_ALIGNMENT", tumorEvidence.avgAltMinDistanceFromAlignment())
                .computeEndFromAlleles(alleles, (int) tumorEvidence.position())
                .source("SAGE")
                .genotypes(tumor, normal)
                .alleles(alleles);

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(tumorEvidence.altQuality() / -10d);
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
        header.addMetaDataLine(new VCFInfoHeaderLine("MAP_Q", UNBOUNDED, VCFHeaderLineType.Float, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("MAP_BASE_Q", UNBOUNDED, VCFHeaderLineType.Float, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("AVG_DISTANCE_RECORD", UNBOUNDED, VCFHeaderLineType.Float, "TODO"));
        header.addMetaDataLine(new VCFInfoHeaderLine("AVG_DISTANCE_ALIGNMENT", UNBOUNDED, VCFHeaderLineType.Float, "TODO"));

        return header;
    }

    @Override
    public void close() {
        writer.close();
    }

}
