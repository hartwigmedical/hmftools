package com.hartwig.hmftools.common.hotspot;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class HotspotEvidenceVCF {


    private final String tumorSample;
    private final String normalSample;

    public HotspotEvidenceVCF(@NotNull final String normalSample, @NotNull final String tumorSample) {
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;
    }

    public void write(@NotNull final String filename, @NotNull final List<HotspotEvidence> evidence) {
        final VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(filename)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .build();

        VCFHeader header = new VCFHeader(Collections.emptySet(), Lists.newArrayList(normalSample, tumorSample));
        header.addMetaDataLine(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read Depth"));
        header.addMetaDataLine(new VCFFormatHeaderLine("AD", 1, VCFHeaderLineType.Integer, "Allelic Depth"));
        header.addMetaDataLine(new VCFInfoHeaderLine("HT", 1, VCFHeaderLineType.String, "Hotspot Type: INFRAME, SNV, MNV, INSERT or DELETE"));
        writer.setHeader(header);
        writer.writeHeader(header);

        for (final HotspotEvidence hotspotEvidence : evidence) {

            final Allele ref = Allele.create(hotspotEvidence.ref(), true);
            final Allele alt = Allele.create(hotspotEvidence.alt(), false);
            final List<Allele> alleles = Lists.newArrayList(ref, alt);

            Genotype tumor = new GenotypeBuilder(tumorSample).DP(hotspotEvidence.tumorReads())
                    .AD(new int[] { hotspotEvidence.tumorEvidence() })
                    .alleles(alleles)
                    .make();

            Genotype normal = new GenotypeBuilder(normalSample).DP(hotspotEvidence.normalReads())
                    .AD(new int[] { hotspotEvidence.normalEvidence() })
                    .alleles(alleles)
                    .make();

            VariantContext context = new VariantContextBuilder().chr(hotspotEvidence.chromosome())
                    .start(hotspotEvidence.position())
                    .attribute("HT", hotspotEvidence.type().toString())
                    .computeEndFromAlleles(tumor.getAlleles(), (int) hotspotEvidence.position())
                    .source(hotspotEvidence.type().toString())
                    .genotypes(tumor, normal)
                    .alleles(alleles)
                    .make();

            context.getCommonInfo().setLog10PError(hotspotEvidence.qualityScore() / -10d);
            writer.add(context);
        }
        writer.close();
    }

}
