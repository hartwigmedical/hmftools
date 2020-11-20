package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.List;

import org.apache.commons.compress.utils.Lists;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

class GermlineWhitelist {

    static List<VariantContext> grch37Whitelist() {
        List<VariantContext> result = Lists.newArrayList();
        result.add(create("3", 10191605, "C", "T"));
        result.add(create("13", 32953884, "C", "T"));
        return result;
    }

    static List<VariantContext> grch38Whitelist() {
        List<VariantContext> result = Lists.newArrayList();
        result.add(create("chr3", 10149921, "C", "T"));
        result.add(create("chr13", 32379747, "C", "T"));
        return result;
    }

    private static VariantContext create(String contig, int position, String ref, String alt) {
        return create(contig, position, ref, alt, GermlineHotspotVCF.WHITELIST_FLAG);
    }

    static VariantContext create(String contig, int position, String ref, String alt, String flag) {
        List<Allele> alleles = Lists.newArrayList();
        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(alt, false));

        return new VariantContextBuilder().alleles(alleles)
                .chr(contig)
                .start(position)
                .computeEndFromAlleles(alleles, position)
                .attribute(flag, true)
                .make();
    }

}
