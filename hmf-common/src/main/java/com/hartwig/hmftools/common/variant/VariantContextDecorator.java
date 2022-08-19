package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.MAPPABILITY_TAG;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.localPhaseSetsStr;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.PURPLE_VARIANT_CN_INFO;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_COUNT_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_SEQUENCE_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.TRINUCLEOTIDE_FLAG;
import static com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser.VAR_IMPACT;

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

    @NotNull
    public PathogenicSummary clinvarPathogenicSummary()
    {
        if(mClinvarPathogenicSummary == null)
        {
            mClinvarPathogenicSummary = PathogenicSummaryFactory.fromContext(mContext);
        }

        return mClinvarPathogenicSummary;
    }

    @NotNull
    public VariantContext context()
    {
        return mContext;
    }

    @NotNull
    public String filter()
    {
        return mFilter;
    }

    @NotNull
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

    @NotNull
    public VariantType type()
    {
        return mType;
    }

    @NotNull
    public String ref()
    {
        return mRef;
    }

    @NotNull
    public String alt()
    {
        return mAlt;
    }

    @NotNull
    public VariantImpact variantImpact()
    {
        if(mVariantImpact == null)
        {
            mVariantImpact = VariantImpactSerialiser.fromVariantContext(mContext);
        }

        return mVariantImpact;
    }

    @NotNull
    public String gene()
    {
        return variantImpact().CanonicalGeneName;
    }

    @NotNull
    public DriverImpact impact()
    {
        if(mDriverImpact == null)
        {
            mDriverImpact = DriverImpact.select(mType, variantImpact().CanonicalCodingEffect);
        }

        return mDriverImpact;
    }

    @NotNull
    public CodingEffect canonicalCodingEffect()
    {
        return variantImpact().CanonicalCodingEffect;
    }

    public double qual()
    {
        return mContext.getPhredScaledQual();
    }

    public double adjustedCopyNumber() { return mContext.getAttributeAsDouble(PURPLE_CN_INFO, 0); }

    public double adjustedVaf()
    {
        return mContext.getAttributeAsDouble(PURPLE_AF_INFO, 0);
    }

    public boolean biallelic() { return mContext.getAttributeAsBoolean(PURPLE_BIALLELIC_FLAG, false); }

    public double minorAlleleCopyNumber() { return mContext.getAttributeAsDouble(PURPLE_MINOR_ALLELE_CN_INFO, 0); }

    public double variantCopyNumber() { return mContext.getAttributeAsDouble(PURPLE_VARIANT_CN_INFO, 0); }

    @Nullable
    public String localPhaseSetsToString()
    {
        List<Integer> localPhaseSets = mContext.getAttributeAsIntList(VariantVcfTags.LOCAL_PHASE_SET, 0);
        return localPhaseSetsStr(localPhaseSets);
    }

    @Nullable
    public Integer localPhaseSet()
    {
        if(!mContext.hasAttribute(VariantVcfTags.LOCAL_PHASE_SET))
            return null;

        List<Integer> localPhaseSets = mContext.getAttributeAsIntList(VariantVcfTags.LOCAL_PHASE_SET, 0);
        return !localPhaseSets.isEmpty() ? localPhaseSets.get(0) : null;
    }

    @NotNull
    public AllelicDepth allelicDepth(@NotNull final String sample)
    {
        final Genotype genotype = mContext.getGenotype(sample);
        return genotype != null ? AllelicDepth.fromGenotype(genotype) : NO_DEPTH;
    }

    @NotNull
    public GenotypeStatus genotypeStatus(@NotNull final String sample)
    {
        final Genotype genotype = mContext.getGenotype(sample);
        return genotype != null ? GenotypeStatus.fromGenotype(genotype) : GenotypeStatus.UNKNOWN;
    }

    @NotNull
    public VariantTier tier()
    {
        return mTier;
    }

    @NotNull
    public PathogenicSummary pathogenicSummary()
    {
        return PathogenicSummaryFactory.fromContext(mContext);
    }

    public int repeatCount()
    {
        return mContext.getAttributeAsInt(REPEAT_COUNT_FLAG, 0);
    }

    @NotNull
    public String repeatSequence()
    {
        return mContext.getAttributeAsString(REPEAT_SEQUENCE_FLAG, Strings.EMPTY);
    }

    @NotNull
    public Hotspot hotspot()
    {
        return Hotspot.fromVariant(mContext);
    }

    public boolean isHotspot()
    {
        return hotspot() == Hotspot.HOTSPOT;
    }

    @NotNull
    public String trinucleotideContext()
    {
        return mContext.getAttributeAsString(TRINUCLEOTIDE_FLAG, Strings.EMPTY);
    }

    public double mappability()
    {
        return mContext.getAttributeAsDouble(MAPPABILITY_TAG, 0);
    }

    public boolean reported()
    {
        return mContext.getAttributeAsBoolean(REPORTED_FLAG, false);
    }

    @NotNull
    public String microhomology()
    {
        return mContext.getAttributeAsString(MICROHOMOLOGY_FLAG, Strings.EMPTY);
    }

    public boolean isPathogenic()
    {
        if(clinvarPathogenicSummary().pathogenicity() == Pathogenicity.BENIGN_BLACKLIST)
        {
            return false;
        }

        if(isHotspot() || clinvarPathogenicSummary().pathogenicity().isPathogenic())
        {
            return true;
        }

        return clinvarPathogenicSummary().pathogenicity() == Pathogenicity.UNKNOWN
                && PATHOGENIC_EFFECT.contains(variantImpact().CanonicalCodingEffect);
    }

    @NotNull
    private static String displayFilter(@NotNull final VariantContext context)
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

    private static final AllelicDepth NO_DEPTH = new AllelicDepth()
    {
        @Override
        public int totalReadCount()
        {
            return 0;
        }

        @Override
        public int alleleReadCount()
        {
            return 0;
        }
    };

    @Override
    public String toString()
    {
        return chromosome() + ":" + position() + " " + mRef + '>' + mAlt;
    }
}
