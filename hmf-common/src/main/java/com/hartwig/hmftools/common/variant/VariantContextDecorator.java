package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.AllelicDepth.NO_DEPTH;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.MAPPABILITY_TAG;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.localPhaseSetsStr;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory;
import com.hartwig.hmftools.common.pathogenic.Pathogenicity;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextDecorator implements GenomePosition
{
    private static final Set<CodingEffect> PATHOGENIC_EFFECT = EnumSet.of(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.SPLICE);

    private final VariantContext mContext;
    private final VariantType mType;
    private final String mFilter;
    private final String mRef;
    private final String mAlt;
    private final VariantTier mTier;

    @Nullable
    private VariantImpact mVariantImpact;
    @Nullable
    private DriverImpact mDriverImpact;
    @Nullable
    private PathogenicSummary mClinvarPathogenicSummary;

    public VariantContextDecorator(final VariantContext context)
    {
        mContext = context;
        mFilter = displayFilter(context);
        mType = VariantType.type(context);
        mRef = getRef(context);
        mAlt = getAlt(context);
        mTier = VariantTier.fromContext(context);
        mVariantImpact = null;
        mDriverImpact = null;
        mClinvarPathogenicSummary = null;
    }

    public static String getRef(final VariantContext context)
    {
        return context.getReference().getBaseString();
    }

    public static String getAlt(final VariantContext context)
    {
        return context.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.joining(","));
    }

    public boolean isPass()
    {
        return mFilter.equals(SomaticVariantFactory.PASS_FILTER);
    }

    public PathogenicSummary clinvarPathogenicSummary()
    {
        if(mClinvarPathogenicSummary == null)
        {
            mClinvarPathogenicSummary = PathogenicSummaryFactory.fromContext(mContext);
        }

        return mClinvarPathogenicSummary;
    }

    public VariantContext context()
    {
        return mContext;
    }

    public String filter()
    {
        return mFilter;
    }

    public VariantType type()
    {
        return mType;
    }

    public String ref()
    {
        return mRef;
    }

    public String alt()
    {
        return mAlt;
    }

    @Override
    public String chromosome()
    {
        return mContext.getContig();
    }

    @Override
    public int position()
    {
        return mContext.getStart();
    }

    public VariantImpact variantImpact()
    {
        if(mVariantImpact == null)
        {
            mVariantImpact = VariantImpactSerialiser.fromVariantContext(mContext);
        }

        return mVariantImpact;
    }

    public String gene()
    {
        return variantImpact().GeneName;
    }

    public DriverImpact impact()
    {
        if(mDriverImpact == null)
        {
            mDriverImpact = DriverImpact.select(mType, variantImpact().CanonicalCodingEffect);
        }

        return mDriverImpact;
    }

    public CodingEffect canonicalCodingEffect()
    {
        return variantImpact().CanonicalCodingEffect;
    }

    public double qual()
    {
        return mContext.getPhredScaledQual();
    }

    public double adjustedCopyNumber()
    {
        return mContext.getAttributeAsDouble(PURPLE_CN, 0);
    }

    public double adjustedVaf()
    {
        return mContext.getAttributeAsDouble(PURPLE_AF, 0);
    }

    public boolean biallelic()
    {
        return mContext.getAttributeAsBoolean(PURPLE_BIALLELIC_FLAG, false);
    }

    public double minorAlleleCopyNumber()
    {
        return mContext.getAttributeAsDouble(PURPLE_MINOR_ALLELE_CN_INFO, 0);
    }

    public double variantCopyNumber()
    {
        return mContext.getAttributeAsDouble(PURPLE_VARIANT_CN, 0);
    }

    @Nullable
    public String localPhaseSetsToString()
    {
        List<Integer> localPhaseSets = mContext.getAttributeAsIntList(LOCAL_PHASE_SET, 0);
        return localPhaseSetsStr(localPhaseSets);
    }

    @Nullable
    public Integer localPhaseSet()
    {
        if(!mContext.hasAttribute(LOCAL_PHASE_SET))
        {
            return null;
        }

        List<Integer> localPhaseSets = mContext.getAttributeAsIntList(LOCAL_PHASE_SET, 0);
        return !localPhaseSets.isEmpty() ? localPhaseSets.get(0) : null;
    }

    public AllelicDepth allelicDepth(final String sample)
    {
        final Genotype genotype = mContext.getGenotype(sample);
        return genotype != null ? AllelicDepth.fromGenotype(genotype) : NO_DEPTH;
    }

    public GenotypeStatus genotypeStatus(final String sample)
    {
        final Genotype genotype = mContext.getGenotype(sample);
        return genotype != null ? GenotypeStatus.fromGenotype(genotype) : GenotypeStatus.UNKNOWN;
    }

    public VariantTier tier()
    {
        return mTier;
    }

    public PathogenicSummary pathogenicSummary()
    {
        return PathogenicSummaryFactory.fromContext(mContext);
    }

    public int repeatCount()
    {
        return mContext.getAttributeAsInt(SageVcfTags.REPEAT_COUNT, 0);
    }

    public String repeatSequence()
    {
        return mContext.getAttributeAsString(SageVcfTags.REPEAT_SEQUENCE, Strings.EMPTY);
    }

    public Hotspot hotspot()
    {
        return Hotspot.fromVariant(mContext);
    }

    public boolean isHotspot()
    {
        return hotspot() == Hotspot.HOTSPOT;
    }

    public String trinucleotideContext()
    {
        return mContext.getAttributeAsString(TRINUCLEOTIDE_CONTEXT, Strings.EMPTY);
    }

    public double mappability()
    {
        return mContext.getAttributeAsDouble(MAPPABILITY_TAG, 0);
    }

    public boolean reported()
    {
        return mContext.getAttributeAsBoolean(REPORTED_FLAG, false);
    }

    public String microhomology()
    {
        return mContext.getAttributeAsString(MICROHOMOLOGY, Strings.EMPTY);
    }

    public boolean isPathogenic()
    {
        if(clinvarPathogenicSummary().Status == Pathogenicity.BENIGN_BLACKLIST)
        {
            return false;
        }

        if(isHotspot() || clinvarPathogenicSummary().Status.isPathogenic())
        {
            return true;
        }

        return clinvarPathogenicSummary().Status == Pathogenicity.UNKNOWN
                && PATHOGENIC_EFFECT.contains(variantImpact().CanonicalCodingEffect);
    }

    private static String displayFilter(final VariantContext context)
    {
        if(context.isFiltered())
        {
            StringJoiner joiner = new StringJoiner(";");
            context.getFilters().forEach(joiner::add);
            return joiner.toString();
        }
        else
        {
            return SomaticVariantFactory.PASS_FILTER;
        }
    }

    @Override
    public String toString()
    {
        return chromosome() + ":" + position() + " " + mRef + '>' + mAlt;
    }
}
