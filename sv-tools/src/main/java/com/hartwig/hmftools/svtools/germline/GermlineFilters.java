package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.IMPRECISE;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class GermlineFilters implements VariantContextFilter
{
    private GermlineVcfConfig mConfig;

    // VCF field identifiers
    public final static String QUAL = "QUAL";
    public final static String SR = "SR";
    public final static String SRQ = "SRQ";
    public final static String VF = "VF";
    public final static String RP = "RP";
    public final static String RPQ = "RPQ";
    public final static String REF = "REF";
    public final static String BEID = "BEID";
    public final static String BEIDL = "BEIDL";
    public final static String HOMSEQ = "HOMSEQ";

    public final static String AS = "AS";
    public final static String CAS = "CAS";
    public final static String RAS = "RAS";

    private final static String BASRP = "BASRP";
    private final static String ASRP = "ASRP";
    private final static String SB = "SB";
    private final static String BVF = "BVF";
    private final static String REFPAIR = "REFPAIR";

    private final static int SHORT_DEL_DUP_LENGTH = 1000;
    private static final String LONG_POLY_C = "CCCCCCCCCCCCCCCC";
    private static final String LONG_POLY_G = "GGGGGGGGGGGGGGGG";
    private static final double MAX_ALLOWABLE_SHORT_EVENT_STRAND_BIAS = 0.95;
    private static final double MIN_AF = 0.005;

    public GermlineFilters(final GermlineVcfConfig config)
    {
        mConfig = config;
    }

    @Override
    public boolean test(final VariantContext variant)
    {
        if (isImprecise(variant))
            return false;

        if (invalidAF(variant))
            return false;

        if (noASRP(variant))
            return false;

        return true;
    }

    public static boolean invalidAF(final VariantContext variant)
    {
        // Filter variants with an allelic fraction of less than 10% in the germline
        final Genotype normalData = variant.getGenotype(0);
        double vf = getDoubleValue(normalData, VF);
        double bvf = getDoubleValue(normalData, BVF);
        double ref = getDoubleValue(normalData, REF);
        double refPair = getDoubleValue(normalData, REFPAIR);
        double af = (vf + bvf) / (vf + bvf + ref + refPair);

        return af < MIN_AF;
    }

    public boolean belowQualScoreThreshold(final VariantContext variant)
    {
        double qualScore = getDoubleValue(variant.getGenotype(0), QUAL);
        return qualScore < mConfig.QualScoreThreshold;
    }

    public static boolean isImprecise(final VariantContext variant)
    {
        return variant.getCommonInfo().getAttribute(IMPRECISE) != null;
    }

    public static boolean noASRP(final VariantContext variant)
    {
        // Filter single breakend variants without an assembly containing at least one discordant read pair
        return variant.getCommonInfo().getAttributeAsInt(BASRP, 0) == 0;
    }

    public static boolean longPolyCorG(final StructuralVariant sv, final VariantContext variant)
    {
        final String insertSeq = sv.insertSequence();
        return insertSeq.contains(LONG_POLY_C) || insertSeq.contains(LONG_POLY_G);
    }

    private static boolean isShortDelDup(final StructuralVariant sv)
    {
        return ((sv.type() == StructuralVariantType.DEL || sv.type() == StructuralVariantType.DUP)
            && (sv.position(false) - sv.position(true) <= SHORT_DEL_DUP_LENGTH));
    }

    public static boolean zeroDiscordantReadSupport(final StructuralVariant sv, final VariantContext variant)
    {
        // Filter breakpoints with no discordant read pair support (either directly, or via assembly)
        // either directly, or through assembly which are not deletion or duplications under 1000bp
        if(!isShortDelDup(sv))
            return false;

        return (variant.getCommonInfo().getAttributeAsInt(ASRP, 0) == 0
        && variant.getCommonInfo().getAttributeAsInt(RP, 0) == 0);
    }

    public static boolean zeroSplitReadSupport(final StructuralVariant sv, final VariantContext variant)
    {
        // Filter breakpoints with no split read support either directly, or through assembly which are not deletion or duplications under 1000bp
        if(!isShortDelDup(sv))
            return false;

        return variant.getCommonInfo().getAttributeAsInt(SR, 0) == 0;
    }

    public static boolean hasStrandBias(final StructuralVariant sv, final VariantContext variant)
    {
        // Filter deletion or duplication breakpoints under 1000bp with a split read strand bias of 0.95 or greater
        if(!isShortDelDup(sv))
            return false;

        return variant.getCommonInfo().getAttributeAsDouble(SB, 0.0) >= MAX_ALLOWABLE_SHORT_EVENT_STRAND_BIAS;
    }

    public static double getDoubleValue(final Genotype data, final String field)
    {
        Object value = data.getAnyAttribute(field);
        return value != null ? Double.parseDouble((String)value) : 0.0;
    }

    public static int getIntValue(final Genotype data, final String field)
    {
        Object value = data.getAnyAttribute(field);
        return value != null ? Integer.parseInt((String)value) : 0;
    }

}
