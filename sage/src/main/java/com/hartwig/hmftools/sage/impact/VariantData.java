package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.variant.VariantConsequence.NON_CODING_TRANSCRIPT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.consequencesToString;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.VariantType.MNP;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.DELIM;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.ITEM_DELIM;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    // other key data
    public boolean mPhasedInframeIndel;
    public String mMicrohomology;
    public String mRepeatSequence;

    private final int mIndelBaseDiff;

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

        mGeneImpacts = Maps.newHashMap();
    }

    public static VariantData fromContext(final VariantContext variantContext)
    {
        int variantPosition = (int)variantContext.getStart();
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

        if(geneImpacts.stream().filter(x -> x.TransData != null) .anyMatch(x -> x.TransData.TransId == impact.TransData.TransId))
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
