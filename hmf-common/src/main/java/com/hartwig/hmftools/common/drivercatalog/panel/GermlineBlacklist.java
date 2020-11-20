package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.List;

import org.apache.commons.compress.utils.Lists;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineBlacklist {


    static List<VariantContext> grch37Blacklist() {
        List<VariantContext> result = Lists.newArrayList();
        // FANCL
        result.add(create("2", 58386928, "G", "GTAAT"));

        // FGFR4
        result.add(create("5", 176520243, "G", "A"));

        // WT1
        result.add(create("11", 32413549, "CT", "C"));
        result.add(create("11", 32413552, "CA", "C"));
        result.add(create("11", 32413561, "GGA", "G"));
        result.add(create("11", 32413573, "C", "CA"));

        // MEN1
        result.add(create("11", 64573217, "C", "CG"));
        result.add(create("11", 64573235, "AGTTG", "A"));

        // TSC2
        result.add(create("16", 2138242, "G", "GTTTTT"));
        result.add(create("16", 2138244, "A", "AG"));
        result.add(create("16", 2138249, "A", "AC"));
        result.add(create("16", 2138253, "GCTCC", "G"));

        return result;
    }

    static VariantContext create(String contig, int position, String ref, String alt) {
        return GermlineWhitelist.create(contig, position, ref, alt, GermlineBlacklistVCF.BLACKLIST_FLAG);

    }

}
