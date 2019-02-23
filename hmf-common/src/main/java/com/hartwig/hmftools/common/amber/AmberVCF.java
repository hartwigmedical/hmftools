package com.hartwig.hmftools.common.amber;

import java.util.Collection;
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
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class AmberVCF {

    private final static String PASS = "PASS";

    private final String tumorSample;
    private final String normalSample;
    private final VCFHeader header;

    public AmberVCF(@NotNull final String normalSample, @NotNull final String tumorSample) {
        this.tumorSample = tumorSample;
        this.normalSample = normalSample;

        this.header = new VCFHeader(Collections.emptySet(), Lists.newArrayList(normalSample, tumorSample));
        header.addMetaDataLine(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Read Depth"));
        header.addMetaDataLine(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Allelic Depth"));
        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));
    }

    public void write(@NotNull final String filename, @NotNull final Collection<TumorBAF> evidence) {
        final List<TumorBAF> list = Lists.newArrayList(evidence);
        Collections.sort(list);

        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        writer.setHeader(header);
        writer.writeHeader(header);

        list.forEach(x -> writer.add(create(x)));
        writer.close();
    }


    public void writeContamination(@NotNull final String filename, @NotNull final Collection<TumorContamination> evidence) {
        final List<TumorContamination> list = Lists.newArrayList(evidence);
        Collections.sort(list);

        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        writer.setHeader(header);
        writer.writeHeader(header);

        list.forEach(x -> writer.add(create(x)));
        writer.close();
    }

    @NotNull
    private VariantContext create(@NotNull final TumorBAF tumorBaf) {

        final Allele ref = Allele.create(tumorBaf.ref(), true);
        final Allele alt = Allele.create(tumorBaf.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final Genotype tumor = new GenotypeBuilder(tumorSample).DP(tumorBaf.tumorReadDepth())
                .AD(new int[] { tumorBaf.tumorRefSupport(), tumorBaf.tumorAltSupport() })
                .alleles(alleles)
                .make();

        final Genotype normal = new GenotypeBuilder(normalSample).DP(tumorBaf.normalReadDepth())
                .AD(new int[] { tumorBaf.normalRefSupport(), tumorBaf.normalAltSupport() })
                .alleles(alleles)
                .make();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(tumorBaf.chromosome())
                .start(tumorBaf.position())
                .computeEndFromAlleles(alleles, (int) tumorBaf.position())
                .genotypes(tumor, normal)
                .alleles(alleles);

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(tumorBaf.tumorAltQuality() / -10d);
        return context;
    }

    @NotNull
    private VariantContext create(@NotNull final TumorContamination contamination) {
        assert(contamination.normal().altSupport() == 0);

        final Allele ref = Allele.create(contamination.tumor().ref().toString(), true);
        final Allele alt = Allele.create(contamination.tumor().alt().toString(), false);

        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final Genotype tumor = new GenotypeBuilder(tumorSample).DP(contamination.tumor().readDepth())
                .AD(new int[] { contamination.tumor().refSupport(), contamination.tumor().altSupport() })
                .alleles(alleles)
                .make();

        final Genotype normal = new GenotypeBuilder(normalSample).DP(contamination.normal().readDepth())
                .AD(new int[] { contamination.normal().refSupport(), contamination.normal().altSupport()})
                .alleles(alleles)
                .make();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(contamination.chromosome())
                .start(contamination.position())
                .computeEndFromAlleles(alleles, (int) contamination.position())
                .genotypes(tumor, normal)
                .alleles(alleles);

        return builder.make();
    }
}
