package com.hartwig.hmftools.amber;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.amber.contamination.TumorContamination;

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
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class VCFWriter
{
    private static final String PASS = "PASS";

    private final AmberConfig mConfig;

    public VCFWriter(final AmberConfig config)
    {
        mConfig = config;
    }

    void writeContamination(final String filename, final Collection<TumorContamination> evidence)
    {
        final List<TumorContamination> list = Lists.newArrayList(evidence);
        Collections.sort(list);

        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, true).build();

        final VCFHeader header = header(Lists.newArrayList(mConfig.primaryReference(), mConfig.TumorId));
        writer.setHeader(header);
        writer.writeHeader(header);

        list.forEach(x -> writer.add(create(x)));
        writer.close();
    }

    public static void writeBaseDepths(final String filename, final Collection<PositionEvidence> baseDepths, String sampleName)
    {
        final List<PositionEvidence> list = Lists.newArrayList(baseDepths);
        Collections.sort(list);

        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        final VCFHeader header = header(Collections.singletonList(sampleName));
        writer.setHeader(header);
        writer.writeHeader(header);

        list.forEach(x -> writer.add(create(x, sampleName)));
        writer.close();
    }

    @NotNull
    private VariantContext create(final TumorContamination contamination)
    {
        final Allele ref = Allele.create(contamination.Tumor.ref().toString(), true);
        final Allele alt = Allele.create(contamination.Tumor.alt().toString(), false);

        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final Genotype tumor = new GenotypeBuilder(mConfig.TumorId).DP(contamination.Tumor.readDepth())
                .AD(new int[] { contamination.Tumor.refSupport(), contamination.Tumor.altSupport() })
                .alleles(alleles)
                .make();

        //        final Genotype normal = new GenotypeBuilder(mConfig.primaryReference()).DP(contamination.Normal.readDepth())
        //                .AD(new int[] { contamination.Normal.refSupport(), contamination.Normal.altSupport() })
        //                .alleles(alleles)
        //                .make();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(contamination.chromosome())
                .start(contamination.position())
                .computeEndFromAlleles(alleles, contamination.position())
                .genotypes(tumor)
                .alleles(alleles);

        return builder.make();
    }

    @NotNull
    private static VariantContext create(final PositionEvidence snp, String sampleName)
    {
        final List<Allele> alleles = Lists.newArrayList();
        alleles.add(Allele.create(snp.ref().toString(), true));
        alleles.add(Allele.create(snp.alt().toString(), false));

        final List<Integer> adField = Lists.newArrayList();
        adField.add(snp.RefSupport);
        adField.add(snp.AltSupport);

        final Genotype normal = new GenotypeBuilder(sampleName)
                .DP(snp.ReadDepth)
                .AD(adField.stream().mapToInt(i -> i).toArray())
                .alleles(alleles)
                .make();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(snp.chromosome())
                .start(snp.position())
                .computeEndFromAlleles(alleles, snp.position())
                .genotypes(normal)
                .alleles(alleles);

        return builder.make();
    }

    private static VCFHeader header(final List<String> samples)
    {
        VCFHeader header = new VCFHeader(Collections.emptySet(), Lists.newArrayList(samples));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));
        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));

        return header;
    }

    private static List<Allele> alleles(final PositionEvidence baseDepth)
    {
        final Allele ref = Allele.create(baseDepth.ref(), true);
        final Allele alt = Allele.create(baseDepth.alt(), false);

        return Lists.newArrayList(ref, alt);
    }

    private static List<Allele> alleles(final TumorBAF baf)
    {
        final Allele ref = Allele.create(baf.ref(), true);
        final Allele alt = Allele.create(baf.alt(), false);

        return Lists.newArrayList(ref, alt);
    }
}
