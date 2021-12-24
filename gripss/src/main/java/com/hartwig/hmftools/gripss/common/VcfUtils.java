package com.hartwig.hmftools.gripss.common;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.CIPOS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.HOTSPOT;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.IMPRECISE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.LOCAL_LINKED_BY;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.MATE_ID;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PAR_ID;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REMOTE_LINKED_BY;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.TAF;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.LINE_POLY_AT_TEST_LEN;

import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.Interval;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VcfUtils
{
    // VCF fields used by Gripss

    // qual-related
    public static final String VT_QUAL = "QUAL";
    public static final String VT_BQ = "BQ";
    public static final String VT_BAQ = "BAQ";
    public static final String VT_SRQ = "SRQ";
    public static final String VT_RPQ = "RPQ";
    public static final String VT_BUMQ = "BUMQ";

    // read counts
    public static final String VT_SR = "SR";
    public static final String VT_VF = "VF";
    public static final String VT_BVF = "BVF";
    public static final String VT_RP = "RP";
    public static final String VT_IC = "IC";
    public static final String VT_REF = "REF";
    public static final String VT_BUM = "BUM";
    public static final String VT_ASRP = "ASRP";
    public static final String VT_ASSR = "ASSR";
    public static final String VT_BASRP = "BASRP";
    public static final String VT_BASSR = "BASSR";

    // other links and info
    public static final String VT_BEID = "BEID";
    public static final String VT_BEIDL = "BEIDL";
    public static final String VT_HOMSEQ = "HOMSEQ";
    public static final String VT_IHOMPOS = "IHOMPOS";

    public static final String VT_PAR_ID = PAR_ID;

    public static final String VT_AS = "AS";
    public static final String VT_CAS = "CAS";
    public static final String VT_RAS = "RAS";

    public static final String VT_EVENT = "EVENT";
    public static final String VT_SB = "SB";
    public static final String VT_BSC = "BSC";
    public static final String VT_REFPAIR = "REFPAIR";
    public static final String VT_CIPOS = CIPOS;
    public static final String VT_CIRPOS = "CIRPOS";
    public static final String VT_IMPRECISE = IMPRECISE;

    public static final String VT_LOCAL_LINKED_BY = LOCAL_LINKED_BY;
    public static final String VT_REMOTE_LINKED_BY = REMOTE_LINKED_BY;
    public static final String VT_TAF = TAF;
    public static final String VT_HOTSPOT = HOTSPOT;
    public static final String VT_REALIGN = "REALIGN";
    public static final String VT_EVENT_TYPE = "EVENTTYPE";
    public static final String VT_ALT_PATH = "ALT_PATH";
    public static final String VT_RESCUE_INFO = "RESCUED";

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
            GR_LOGGER.error("missing sample names in VCF: {}", vcfSampleNames);
            return null;
        }


        return new GenotypeIds(referenceOrdinal, tumorOrdinal, vcfRefefenceId, vcfTumorId);
    }

    public static int getGenotypeAttributeAsInt(final Genotype genotype, final String attribute, int defaultVaue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultVaue : Integer.parseInt(value.toString());
    }

    public static double getGenotypeAttributeAsDouble(final Genotype genotype, final String attribute, double defaultVaue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultVaue : Double.parseDouble(value.toString());
    }

    public static final Interval confidenceInterval(final VariantContext variantContext, final String attribute)
    {
        if(!variantContext.hasAttribute(attribute))
            return new Interval();

        List<Integer> values = variantContext.getAttributeAsIntList(attribute, 0);
        return new Interval(values.get(0), values.get(1));
    }

    public static int sglFragmentCount(final Genotype genotype)
    {
        int bsc = getGenotypeAttributeAsInt(genotype, VT_BSC, 0);
        int basrp = getGenotypeAttributeAsInt(genotype, VT_BASRP, 0);
        int bassr = getGenotypeAttributeAsInt(genotype, VT_BASSR, 0);
        int bvf = getGenotypeAttributeAsInt(genotype, VT_BVF, 0);

        if(bsc == 0 && basrp == 0 && bassr == 0)
            return 0;
        else
            return bvf;
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

    public static boolean isMobileLineElement(final byte orientation, final String insertSequence)
    {
        int insSeqLength = insertSequence.length();
        if(insSeqLength < LINE_POLY_AT_REQ)
            return false;

        final char polyATChar = orientation == POS_ORIENT ? 'T' : 'A';

        int testLength = min(LINE_POLY_AT_TEST_LEN, insSeqLength);
        int allowedNonRequiredChars = testLength - LINE_POLY_AT_REQ;

        for(int i = 0; i < testLength; ++i)
        {
            if(orientation == POS_ORIENT)
            {
                if(insertSequence.charAt(i) != polyATChar)
                    --allowedNonRequiredChars;
            }
            else
            {
                if(insertSequence.charAt(insSeqLength - i - 1) != polyATChar)
                    --allowedNonRequiredChars;
            }

            if(allowedNonRequiredChars < 0)
                return false;
        }

        return true;
    }
}
