package com.hartwig.hmftools.pave;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.VariantType.MNP;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;

import static htsjdk.variant.vcf.VCFConstants.ALLELE_FREQUENCY_KEY;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.apache.logging.log4j.util.Strings;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantData
{
    public final String Chromosome;

    public final int Position; // as per the variant definition

    // position and end position should be symmetrical, so that EndPosition functions like Position on the negative strand
    // pos for an SNV
    // pos + alt-length - 1 for MNV
    // pos + ref-length - 1 for DEL
    // pos + 1 for an INS
    public final int EndPosition;

    public final String Ref;

    public final String Alt;

    // other key data
    private int mLocalPhaseSetId;
    private String mMicrohomology;
    private int mRepeatCount;
    private String mRepeatSequence;

    private final int mIndelBaseDiff;

    // compute and cache base positions which have changed due to this variant:
    // for an insert, none
    // for a delete, the deleted positions, so starting after the first (ref) position
    // for SNVs and MNVs, all of them
    private final List<Integer> mAltPositions;

    private final Map<String,List<VariantTransImpact>> mGeneImpacts;

    private VariantData mRealignedVariant;
    private boolean mIsRealignedVariant;

    // associated data
    private VariantContext mVariantContext;

    private boolean mReportable;

    public final Set<String> mFilters;
    public int mPonSampleCount;
    public int mPonMaxReadCount;
    public int mPonMeanReadCount;
    private Double mGnomadFrequency;

    public static final int NO_LOCAL_PHASE_SET = -1;

    public VariantData(final String chromosome, final int position, final String ref, final String alt)
    {
        Chromosome = chromosome;

        Position = position;
        Ref = ref;
        Alt = alt;

        mVariantContext = null;

        mLocalPhaseSetId = -1;
        mMicrohomology = "";
        mRepeatSequence = "";
        mRepeatCount = 0;

        mIndelBaseDiff = Alt.length() - Ref.length();

        if(mIndelBaseDiff != 0)
        {
            // rare case where an INDEL has the ref base change as well, eg CT -> A
            boolean refBaseChanged = alt.charAt(0) != ref.charAt(0);

            if(isInsert())
            {
                if(refBaseChanged)
                {
                    mAltPositions = Lists.newArrayListWithExpectedSize(1);
                    mAltPositions.add(Position);
                }
                else
                {
                    mAltPositions = Lists.newArrayListWithExpectedSize(0);
                }

                EndPosition = Position + 1;
            }
            else
            {
                // a deletion
                int count = abs(mIndelBaseDiff);

                int altPositions = refBaseChanged ? count + 1 : count;
                mAltPositions = Lists.newArrayListWithExpectedSize(altPositions);
                int altStartIndex = refBaseChanged ? 0 : 1;

                for(int i = altStartIndex; i < Ref.length(); ++i)
                {
                    mAltPositions.add(Position + i);
                }

                EndPosition = Position + count + 1; // the first based after the deleted section
            }
        }
        else
        {
            int count = Ref.length();
            mAltPositions = Lists.newArrayListWithExpectedSize(count);

            for(int i = 0; i < Ref.length(); ++i)
            {
                mAltPositions.add(Position + i);
            }

            EndPosition = Position + count - 1;
        }

        mRealignedVariant = null;
        mIsRealignedVariant = false;
        mGeneImpacts = Maps.newHashMap();
        mReportable = false;

        mFilters = Sets.newHashSet();
        mPonSampleCount = 0;
        mPonMaxReadCount = 0;
        mPonMeanReadCount = 0;
        mGnomadFrequency = null;
    }

    public static VariantData fromContext(final VariantContext variantContext)
    {
        int variantPosition = variantContext.getStart();
        String chromosome = variantContext.getContig();

        String ref = variantContext.getReference().getBaseString();

        // only support the first of multiple alts (as can be the case for Strelka)
        String alt = !variantContext.getAlternateAlleles().isEmpty() ? variantContext.getAlternateAlleles().get(0).toString() : ref;

        if(alt.equals("*") || alt.equals("N")) // unhandled for now
            alt = ref;

        VariantData variant = new VariantData(chromosome, variantPosition, ref, alt);
        variant.setContext(variantContext);

        List<Integer> localPhaseSets = variantContext.getAttributeAsIntList(LOCAL_PHASE_SET, NO_LOCAL_PHASE_SET);

        variant.setVariantDetails(
                !localPhaseSets.isEmpty() ? localPhaseSets.get(0) : NO_LOCAL_PHASE_SET,
                variantContext.getAttributeAsString(MICROHOMOLOGY, Strings.EMPTY),
                variantContext.getAttributeAsString(REPEAT_SEQUENCE, Strings.EMPTY),
                variantContext.getAttributeAsInt(REPEAT_COUNT, 0));

        return variant;
    }

    public VariantType type()
    {
        if(mIndelBaseDiff == 0)
            return Ref.length() == 1 ? SNP : MNP;

        return INDEL;
    }

    public int baseDiff() { return mIndelBaseDiff; }
    public boolean isBaseChange() { return mIndelBaseDiff == 0; }
    public boolean isIndel() { return mIndelBaseDiff != 0; }
    public boolean isInsert() { return mIndelBaseDiff > 0; }
    public boolean isDeletion() { return mIndelBaseDiff < 0; }
    public boolean isMixed() { return mIndelBaseDiff != 0 && Alt.charAt(0) != Ref.charAt(0); }
    public boolean isSnv() { return mIndelBaseDiff == 0 && Ref.length() == 1; }
    public boolean isMnv() { return mIndelBaseDiff == 0 && Ref.length() > 1; }

    public List<Integer> altPositions() { return mAltPositions; }

    public VariantData realignedVariant() { return mRealignedVariant; }

    public void setRealignedVariant(final VariantData variant)
    {
        if(variant != null)
        {
            mRealignedVariant = variant;
            variant.markRealignedVariant();
        }
    }

    public boolean isRealignedVariant() { return mIsRealignedVariant; }
    public void markRealignedVariant() { mIsRealignedVariant = true; }

    public VariantContext context() { return mVariantContext; }
    public void setContext(final VariantContext context) { mVariantContext = context; }

    public boolean reported() { return mReportable; }
    public void markReported() { mReportable = true; }

    public void setVariantDetails(int localPhaseSet, final String microHomology, final String repeatSequece, final int repeatCount)
    {
        mLocalPhaseSetId = localPhaseSet;
        mMicrohomology = microHomology;
        mRepeatSequence = repeatSequece;
        mRepeatCount = repeatCount;
    }

    public boolean hasLocalPhaseSet() { return mLocalPhaseSetId != NO_LOCAL_PHASE_SET; }
    public int localPhaseSet() { return mLocalPhaseSetId; }
    public String microhomology() { return mMicrohomology; }
    public String repeatSequence() { return mRepeatSequence; }
    public int repeatCount() { return mRepeatCount; }

    public Map<String,List<VariantTransImpact>> getImpacts() { return mGeneImpacts; }

    public List<VariantTransImpact> getGeneImpacts(final String geneName) { return mGeneImpacts.get(geneName); }

    public VariantTransImpact getCanonicalTransImpacts(final String geneName)
    {
        List<VariantTransImpact> transImpacts = mGeneImpacts.get(geneName);
        return transImpacts != null ? transImpacts.stream().filter(x -> x.TransData.IsCanonical).findFirst().orElse(null) : null;
    }

    public void addImpact(final String geneName, final VariantTransImpact impact)
    {
        if(impact.effects().isEmpty())
            return;

        List<VariantTransImpact> geneImpacts = mGeneImpacts.get(geneName);

        if(geneImpacts == null)
        {
            geneImpacts = Lists.newArrayList(impact);
            mGeneImpacts.put(geneName, geneImpacts);
            return;
        }

        if(geneImpacts.stream().filter(x -> x.TransData != null).anyMatch(x -> x.TransData.TransId == impact.TransData.TransId))
            return;

        geneImpacts.add(impact);
    }

    public boolean altBasesAbove(int position) { return altBasesOutsidePosition(position, true); }
    public boolean altBasesBelow(int position) { return altBasesOutsidePosition(position, false); }

    public boolean altPositionsWithin(int lowerPos, int upperPos)
    {
        if(isInsert())
            return positionsWithin(Position, EndPosition, lowerPos, upperPos);

        // test if all alt positions are at or inside the limits
        return altBasesAbove(lowerPos - 1) && altBasesBelow(upperPos + 1);
    }

    public boolean altPositionsOverlap(int lowerPos, int upperPos)
    {
        // test if all alt positions are at or inside the limits
        return !altBasesBelow(lowerPos) && !altBasesAbove(upperPos);
    }

    public boolean altBasesOutsidePosition(int position, boolean isHigher)
    {
        // rather than just testing variant position vs a genomic / transcript position, it is sometimes required to know if an altered
        // bases is at or beyond a position
        if(mIndelBaseDiff == 0)
        {
            return isHigher ? Position > position : EndPosition < position;
        }

        if(isInsert())
        {
            return isHigher ? EndPosition > position : Position < position;
        }
        else
        {
            // first base of DEL is the ref
            return isHigher ? Position + 1 > position : EndPosition - 1 < position;
        }
    }

    public VariantTransImpact getRealignedImpact(final String geneName, final VariantTransImpact transImpact)
    {
        if(mRealignedVariant == null)
            return null;

        List<VariantTransImpact> raImpacts = mRealignedVariant.getImpacts().get(geneName);

        if(raImpacts == null)
            return null;

        return raImpacts.stream().filter(x -> x.TransData.TransName.equals(transImpact.TransData.TransName)).findFirst().orElse(null);
    }

    public VariantTier tier()
    {
        if(mVariantContext == null)
            return VariantTier.UNKNOWN;

        return VariantTier.fromContext(mVariantContext);
    }

    public double sampleVaf(final String sampleId)
    {
        Genotype sampleGenotype = mVariantContext.getGenotype(sampleId);
        return sampleGenotype != null ? getGenotypeAttributeAsDouble(sampleGenotype, ALLELE_FREQUENCY_KEY, 0) : 0;
    }

    public Set<String> filters() { return mFilters; }
    public void addFilter(final String filter) { mFilters.add(filter); }

    public int ponSampleCount() { return mPonSampleCount; }
    public int ponMaxReadCount() { return mPonMaxReadCount; }
    public int ponMeanReadCount() { return mPonMaxReadCount; }

    public void setPonFrequency(int sampleCount, int maxReadCount, int meanReadCount)
    {
        mPonSampleCount = sampleCount;
        mPonMaxReadCount = maxReadCount;
        mPonMeanReadCount = meanReadCount;
    }

    public Double gnomadFrequency() { return mGnomadFrequency; }
    public void setGnomadFrequency(Double frequency) { mGnomadFrequency = frequency; }

    public String toString()
    {
        if(mIndelBaseDiff == 0 && Ref.length() == 1)
            return String.format("pos(%s:%d) variant(%s>%s)", Chromosome, Position, Ref, Alt);
        else
            return String.format("pos(%s:%d-%d) variant(%s>%s)", Chromosome, Position, EndPosition, Ref, Alt);
    }

    public static String tsvHeader()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("Chromosome");
        sj.add("Position");
        sj.add("Type");
        sj.add("Ref");
        sj.add("Alt");
        sj.add("LPS");
        return sj.toString();
    }

    public String tsvData()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(Chromosome);
        sj.add(String.valueOf(Position));
        sj.add(String.valueOf(type()));
        sj.add(Ref);
        sj.add(Alt);
        sj.add(String.valueOf(mLocalPhaseSetId));

        return sj.toString();
    }

    public static String extraDataHeader()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("Microhomology");
        sj.add("RepeatSequence");
        sj.add("RepeatCount");
        sj.add("LocalPhaseSetId");
        sj.add("Filter");
        sj.add("Tier");
        sj.add("Qual");
        sj.add("AlleleReadCount");
        sj.add("TotalReadCount");
        return sj.toString();
    }

    public String extraDataTsv(final String sampleId)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(mMicrohomology);
        sj.add(mRepeatSequence);
        sj.add(String.valueOf(mRepeatCount));
        sj.add(String.valueOf(mLocalPhaseSetId));

        VariantContextDecorator varDetails = new VariantContextDecorator(mVariantContext);
        AllelicDepth allelicDepth = varDetails.allelicDepth(sampleId);

        sj.add(varDetails.filter());
        sj.add(varDetails.tier().toString());
        sj.add(String.valueOf(varDetails.qual()));
        sj.add(String.valueOf(allelicDepth.AlleleReadCount));
        sj.add(String.valueOf(allelicDepth.TotalReadCount));

        return sj.toString();
    }
}
