package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.svtools.germline.GermlineUtils.GM_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VcfUtils
{
    // VCF fields used by Gripss
    public static final String QUAL = "QUAL";
    public static final String SR = "SR";
    public static final String BQ = "BQ";
    public static final String SRQ = "SRQ";
    public static final String VF = "VF";
    public static final String RP = "RP";
    public static final String IC = "IC";
    public static final String RPQ = "RPQ";
    public static final String REF = "REF";
    public static final String BEID = "BEID";
    public static final String BEIDL = "BEIDL";
    public static final String HOMSEQ = "HOMSEQ";

    public static final String AS = "AS";
    public static final String CAS = "CAS";
    public static final String RAS = "RAS";

    public static final String EVENT = "EVENT";
    public static final String ASRP = "ASRP";
    public static final String SB = "SB";
    public static final String BVF = "BVF";
    public static final String REFPAIR = "REFPAIR";
    public static final String CIRPOS = "CIRPOS";
    public static final String REALIGN = "REALIGN";

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

        int assemblyCount = variantContext.getAttributeAsInt(AS, 0)
                + variantContext.getAttributeAsInt(RAS, 0)
                + variantContext.getAttributeAsInt(CAS, 0);

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
