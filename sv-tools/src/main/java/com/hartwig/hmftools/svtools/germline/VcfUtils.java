package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.CIPOS;
import static com.hartwig.hmftools.svtools.germline.GermlineUtils.GM_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VcfUtils
{
    // VCF fields used by Gripss
    public static final String VT_QUAL = "QUAL";
    public static final String VT_SR = "SR";
    public static final String VT_BQ = "BQ";
    public static final String VT_SRQ = "SRQ";
    public static final String VT_VF = "VF";
    public static final String VT_RP = "RP";
    public static final String VT_IC = "IC";
    public static final String VT_RPQ = "RPQ";
    public static final String VT_REF = "REF";
    public static final String VT_BEID = "BEID";
    public static final String VT_BEIDL = "BEIDL";
    public static final String VT_HOMSEQ = "HOMSEQ";

    public static final String VT_AS = "AS";
    public static final String VT_CAS = "CAS";
    public static final String VT_RAS = "RAS";

    public static final String VT_EVENT = "EVENT";
    public static final String VT_ASRP = "ASRP";
    public static final String VT_SB = "SB";
    public static final String VT_BVF = "BVF";
    public static final String VT_REFPAIR = "REFPAIR";
    public static final String VT_CIPOS = CIPOS;
    public static final String VT_CIRPOS = "CIRPOS";
    public static final String VT_REALIGN = "REALIGN";

    public static GenotypeIds parseVcfSampleIds(final VCFHeader header, final String referenceId, final String tumorId)
    {
        List<String> vcfSampleNames = header.getGenotypeSamples();

        int tumorOrdinal = -1;
        int referenceOrdinal = -1;
        String vcfTumorId = "";
        String vcfRefefenceId = "";

        for(int i = 0; i < vcfSampleNames.size(); ++i)
        {
            String vcfSampleName = vcfSampleNames.get(i);

            if(vcfSampleName.contains(tumorId))
            {
                vcfTumorId = vcfSampleNames.get(i);
                tumorOrdinal = i;
            }
            else if(!referenceId.isEmpty() && vcfSampleName.contains(referenceId))
            {
                vcfRefefenceId = vcfSampleNames.get(i);
                referenceOrdinal = i;
            }
        }

        if(tumorOrdinal < 0 || (!referenceId.isEmpty() && referenceOrdinal < 0))
        {
            GM_LOGGER.error("missing sample names in VCF: {}", vcfSampleNames);
            return null;
        }


        return new GenotypeIds(referenceOrdinal, tumorOrdinal, vcfRefefenceId, vcfTumorId);
    }

    public static int getGenotypeAttributeAsInt(final Genotype genotype, final String attribute, int defaultVaue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultVaue : (int)value;
    }

    public static double getGenotypeAttributeAsDouble(final Genotype genotype, final String attribute, double defaultVaue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultVaue : (double)value;
    }

    public static final Interval confidenceInterval(final VariantContext variantContext, final String attribute)
    {
        if(!variantContext.hasAttribute(attribute))
            return new Interval();

        List<Integer> values = variantContext.getAttributeAsIntList(attribute, 0);
        return new Interval(values.get(0), values.get(1));
    }

    public static List<String> parseAssemblies(final VariantContext variantContext)
    {
        List<String> assemblies = Lists.newArrayList();

        int assemblyCount = variantContext.getAttributeAsInt(VT_AS, 0)
                + variantContext.getAttributeAsInt(VT_RAS, 0)
                + variantContext.getAttributeAsInt(VT_CAS, 0);

        if(assemblyCount >= 2 && variantContext.hasAttribute(VT_BEID) && variantContext.hasAttribute(VT_BEIDL))
        {
            List<String> beids = variantContext.getAttributeAsStringList(VT_BEID, "");
            List<String> beidls = variantContext.getAttributeAsStringList(VT_BEIDL, "");

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
