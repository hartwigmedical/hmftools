package com.hartwig.hmftools.gripss.common;

import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.BEID;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.BEIDL;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_AS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BASRP;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BASSR;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BSC;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_CAS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_RAS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SGL_FRAG_COUNT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VcfUtils
{
    // VCF fields used by Gripss

    public static final Interval confidenceInterval(final VariantContext variantContext, final String attribute)
    {
        if(!variantContext.hasAttribute(attribute))
            return new Interval();

        List<Integer> values = variantContext.getAttributeAsIntList(attribute, 0);
        return new Interval(values.get(0), values.get(1));
    }

    public static int sglFragmentCount(final Genotype genotype)
    {
        int bsc = getGenotypeAttributeAsInt(genotype, GRIDSS_BSC, 0);
        int basrp = getGenotypeAttributeAsInt(genotype, GRIDSS_BASRP, 0);
        int bassr = getGenotypeAttributeAsInt(genotype, GRIDSS_BASSR, 0);
        int bvf = getGenotypeAttributeAsInt(genotype, SGL_FRAG_COUNT, 0);

        if(bsc == 0 && basrp == 0 && bassr == 0)
            return 0;
        else
            return bvf;
    }

    public static List<String> parseAssemblies(final VariantContext variantContext)
    {
        List<String> assemblies = Lists.newArrayList();

        int assemblyCount = variantContext.getAttributeAsInt(GRIDSS_AS, 0)
                + variantContext.getAttributeAsInt(GRIDSS_RAS, 0)
                + variantContext.getAttributeAsInt(GRIDSS_CAS, 0);

        if(assemblyCount >= 2 && variantContext.hasAttribute(BEID) && variantContext.hasAttribute(BEIDL))
        {
            List<String> beids = variantContext.getAttributeAsStringList(BEID, "");
            List<String> beidls = variantContext.getAttributeAsStringList(BEIDL, "");

            if(beidls.size() == beids.size())
            {
                for(int i = 0; i < beids.size(); ++i)
                {
                    assemblies.add(String.format("%s/%s", beids.get(i), beidls.get(i)));
                }
            }
        }

        return assemblies;
    }
}
