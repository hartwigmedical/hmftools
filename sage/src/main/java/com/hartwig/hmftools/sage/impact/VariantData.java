package com.hartwig.hmftools.sage.impact;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.VariantType.MNP;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.DELIM;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.VariantType;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantData
{
    public final String Chromosome;
    public final int Position;
    public final int EndPosition;
    public final String Ref;
    public final String Alt;

    // other key data
    public boolean mPhasedInframeIndel;
    public String mMicrohomology;
    public String mRepeatSequence;

    private final int mIndelBaseDiff;

    // compute and cache base positions which have changed due to this variant:
    // for an insert, none
    // for a delete, the deleted positions, so starting after the first (ref) position
    // for SNVs and MNVs, all of them

    private final List<Integer> mNonRefPositions;

    private final Map<String,List<VariantTransImpact>> mGeneImpacts;

    public VariantData(final String chromosome, final int position, final String ref, final String alt)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;

        mPhasedInframeIndel = false;
        mMicrohomology = "";
        mRepeatSequence = "";

        mIndelBaseDiff = Alt.length() - Ref.length();

        if(isInsert())
        {
            mNonRefPositions = Lists.newArrayListWithExpectedSize(0);
            EndPosition = Position + 1;
        }
        else
        {
            int count = mIndelBaseDiff == 0 ? Ref.length() : abs(mIndelBaseDiff);
            mNonRefPositions = Lists.newArrayListWithExpectedSize(count);

            int startPos = isIndel() ? 1 : 0;
            for(int i = startPos; i < Ref.length(); ++i)
            {
                mNonRefPositions.add(Position + i);
            }

            EndPosition = mNonRefPositions.get(mNonRefPositions.size() - 1);
        }

        mGeneImpacts = Maps.newHashMap();
    }

    public static VariantData fromContext(final VariantContext variantContext)
    {
        int variantPosition = variantContext.getStart();
        String chromosome = variantContext.getContig();

        String ref = variantContext.getReference().getBaseString();
        String alt = variantContext.getAlternateAlleles().stream().map(Allele::toString).collect(Collectors.joining(","));

        return new VariantData(chromosome, variantPosition, ref, alt);
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

    public List<Integer> nonRefPositions() { return mNonRefPositions; }

    public boolean phasedInframeIndel() { return mPhasedInframeIndel; }

    public void setVariantDetails(boolean piIndel, final String mh, final String reqSeq)
    {
        mPhasedInframeIndel = piIndel;
        mMicrohomology = mh;
        mRepeatSequence = reqSeq;
    }

    public String microhomology() { return mMicrohomology; }
    public String repeatSequence() { return mRepeatSequence; }

    public Map<String,List<VariantTransImpact>> getImpacts() { return mGeneImpacts; }

    public void addImpact(final String geneName, final VariantTransImpact impact)
    {
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

    public String toString()
    {
        return String.format("pos(%s:%d) variant(%s>%s)", Chromosome, Position, Ref, Alt);
    }

    public static String csvCommonHeader()
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add("Chromosome").add("Position").add("Type").add("Ref").add("Alt");
        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add(Chromosome);
        sj.add(String.valueOf(Position));
        sj.add(String.valueOf(type()));
        sj.add(Ref);
        sj.add(Alt);

        return sj.toString();
    }


}
