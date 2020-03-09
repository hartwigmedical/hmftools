package com.hartwig.hmftools.amber;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.errorprone.annotations.Var;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSiteFactory;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.TumorContamination;
import com.hartwig.hmftools.common.variant.Variant;

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

public class AmberVCF {

    private final static String PASS = "PASS";

    private final String tumorSample;
    private final String normalSample;
    private final AmberConfig config;

    AmberVCF(@NotNull final AmberConfig config) {
        this.tumorSample = config.tumor();
        this.normalSample = config.reference().get(0);
        this.config = config;
    }

    public void writeBAF(@NotNull final String filename, @NotNull final Collection<TumorBAF> evidence) {
        final List<TumorBAF> list = Lists.newArrayList(evidence);
        Collections.sort(list);

        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        final VCFHeader header = header(config.tumorOnly() ? Collections.singletonList(config.tumor()) : config.allSamples());
        writer.setHeader(header);
        writer.writeHeader(header);

        list.forEach(x -> writer.add(create(x)));
        writer.close();
    }

    public void writeBAF(@NotNull final String filename, @NotNull final Collection<TumorBAF> evidence, List<Collection<BaseDepth>> normal) {
        final List<TumorBAF> list = Lists.newArrayList(evidence);
        Collections.sort(list);

        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        final VCFHeader header = header(config.tumorOnly() ? Collections.singletonList(config.tumor()) : config.allSamples());
        writer.setHeader(header);
        writer.writeHeader(header);

        // Map normals to sites
        final List<Map<AmberSite, BaseDepth>> normalList = Lists.newArrayList();
        for (Collection<BaseDepth> aNormal : normal) {
            Map<AmberSite, BaseDepth> normalMap = Maps.newHashMap();
            for (BaseDepth baseDepth : aNormal) {
                normalMap.put(AmberSiteFactory.asSite(baseDepth), baseDepth);
            }
            normalList.add(normalMap);
        }

        for (TumorBAF tumorBAF : list) {
            AmberSite tumorSite = AmberSiteFactory.tumorSite(tumorBAF);
            List<BaseDepth> normalDepths = Lists.newArrayList();
            for (Map<AmberSite, BaseDepth> normalMap : normalList) {
                normalDepths.add(normalMap.get(tumorSite));
            }

            writer.add(create(tumorBAF, normalDepths));
        }

        writer.close();
    }

    void writeContamination(@NotNull final String filename, @NotNull final Collection<TumorContamination> evidence) {
        final List<TumorContamination> list = Lists.newArrayList(evidence);
        Collections.sort(list);

        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        final VCFHeader header = header(Lists.newArrayList(config.primaryReference(), config.tumor()));
        writer.setHeader(header);
        writer.writeHeader(header);

        list.forEach(x -> writer.add(create(x)));
        writer.close();
    }

    void writeSNPCheck(@NotNull final String filename, @NotNull final List<BaseDepth> baseDepths) {
        final List<BaseDepth> list = Lists.newArrayList(baseDepths);
        Collections.sort(list);

        final VariantContextWriter writer =
                new VariantContextWriterBuilder().setOutputFile(filename).modifyOption(Options.INDEX_ON_THE_FLY, true).build();
        final VCFHeader header = header(Lists.newArrayList(config.primaryReference()));
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
    private VariantContext create(@NotNull final TumorBAF tumorBaf, final List<BaseDepth> normals) {
        final Allele ref = Allele.create(tumorBaf.ref(), true);
        final Allele alt = Allele.create(tumorBaf.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final Genotype tumor = new GenotypeBuilder(tumorSample).DP(tumorBaf.tumorReadDepth())
                .AD(new int[] { tumorBaf.tumorRefSupport(), tumorBaf.tumorAltSupport() })
                .alleles(alleles)
                .make();

        List<Genotype> genotypes = Lists.newArrayList(tumor);
        for (int i = 0; i < normals.size(); i++) {
            final Genotype normal = create(config.reference().get(i), alleles, normals.get(i));
            genotypes.add(normal);
        }

        final VariantContextBuilder builder = new VariantContextBuilder().chr(tumorBaf.chromosome())
                .start(tumorBaf.position())
                .computeEndFromAlleles(alleles, (int) tumorBaf.position())
                .genotypes(genotypes)
                .alleles(alleles);

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(tumorBaf.tumorAltQuality() / -10d);
        return context;
    }


    @NotNull
    private static Genotype create(@NotNull final String sample, @NotNull final List<Allele> alleles, @NotNull final BaseDepth depth) {
        return new GenotypeBuilder(sample).DP(depth.readDepth())
                .AD(new int[] { depth.refSupport(), depth.altSupport() })
                .alleles(alleles)
                .make();
    }

    @NotNull
    private VariantContext create(@NotNull final TumorContamination contamination) {
        assert (contamination.normal().altSupport() == 0);

        final Allele ref = Allele.create(contamination.tumor().ref().toString(), true);
        final Allele alt = Allele.create(contamination.tumor().alt().toString(), false);

        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        final Genotype tumor = new GenotypeBuilder(tumorSample).DP(contamination.tumor().readDepth())
                .AD(new int[] { contamination.tumor().refSupport(), contamination.tumor().altSupport() })
                .alleles(alleles)
                .make();

        final Genotype normal = new GenotypeBuilder(normalSample).DP(contamination.normal().readDepth())
                .AD(new int[] { contamination.normal().refSupport(), contamination.normal().altSupport() })
                .alleles(alleles)
                .make();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(contamination.chromosome())
                .start(contamination.position())
                .computeEndFromAlleles(alleles, (int) contamination.position())
                .genotypes(tumor, normal)
                .alleles(alleles);

        return builder.make();
    }

    @NotNull
    private VariantContext create(@NotNull final BaseDepth snp) {
        final List<Allele> alleles = Lists.newArrayList();
        alleles.add(Allele.create(snp.ref().toString(), true));
        alleles.add(Allele.create(snp.alt().toString(), false));

        final List<Integer> adField = Lists.newArrayList();
        adField.add(snp.refSupport());
        adField.add(snp.altSupport());

        final Genotype normal = new GenotypeBuilder(normalSample).DP(snp.readDepth())
                .AD(adField.stream().mapToInt(i -> i).toArray())
                .alleles(alleles)
                .make();

        final VariantContextBuilder builder = new VariantContextBuilder().chr(snp.chromosome())
                .start(snp.position())
                .computeEndFromAlleles(alleles, (int) snp.position())
                .genotypes(normal)
                .alleles(alleles);

        return builder.make();
    }

    @NotNull
    private static VCFHeader header(final List<String> samples) {
        VCFHeader header = new VCFHeader(Collections.emptySet(), Lists.newArrayList(samples));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));
        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));

        return header;
    }
}
