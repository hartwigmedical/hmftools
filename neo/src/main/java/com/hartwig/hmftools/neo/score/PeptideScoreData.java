package com.hartwig.hmftools.neo.score;

import static java.lang.String.format;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.neo.PeptideData;
import com.hartwig.hmftools.neo.bind.BindData;

public class PeptideScoreData extends PeptideData
{
    private final Set<Integer> mNeIds;

    private double mEffectiveTpm;

    private final List<BindData> mAlleleScoreData;

    private boolean mTpmCalculated;
    private boolean mWritten;

    public PeptideScoreData(final PeptideData peptideData, final int neId)
    {
        super(peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank);

        mNeIds = Sets.newHashSet(neId);

        mEffectiveTpm = 0;

        mAlleleScoreData = Lists.newArrayList();
        mWritten = false;
        mTpmCalculated = false;
    }

    public void addAllele(final String allele)
    {
        BindData bindData = new BindData(allele, Peptide, "", UpFlank, DownFlank);
        bindData.setTPM(effectiveTpm());
        mAlleleScoreData.add(bindData);
    }

    public List<BindData> alleleScoreData() { return mAlleleScoreData; }

    public void addNeoId(int neId) { mNeIds.add(neId); }
    public Set<Integer> neIds() { return mNeIds; }
    public boolean hasNeo(int neId) { return mNeIds.contains(neId); }

    public String neoIdsStr()
    {
        if(mNeIds.size() == 1)
            return String.valueOf(mNeIds.iterator().next());

        StringJoiner sj = new StringJoiner(";");
        mNeIds.forEach(x -> sj.add(String.valueOf(x)));
        return sj.toString();
    }

    public boolean written() { return mWritten; }
    public void setWritten() { mWritten = true; }

    public double effectiveTpm() { return mEffectiveTpm; }
    public boolean tpmCalculated() { return mTpmCalculated; }

    public void setEffectiveTpm(double effectiveTpm)
    {
        mEffectiveTpm = effectiveTpm;
        mTpmCalculated = true;
    }

    public String toString()
    {
        return format("%s neIds(%d) effective(%4.3e)", Peptide, mNeIds.size(), mEffectiveTpm);
    }
}
