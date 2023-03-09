package com.hartwig.hmftools.neo.scorer;

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

    private double mExpectedTpm;
    private double mEffectiveTpm;
    private double mRawEffectiveTpm;

    private final List<BindData> mAlleleScoreData;

    private boolean mTpmCalculated;
    private boolean mWritten;

    public PeptideScoreData(final PeptideData peptideData, final int neId)
    {
        super(peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank);

        mNeIds = Sets.newHashSet(neId);

        mExpectedTpm = 0;
        mEffectiveTpm = 0;
        mRawEffectiveTpm = 0;

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

    public double expectedTpm() { return mExpectedTpm; }
    public double effectiveTpm() { return mEffectiveTpm; }
    public double rawEffectiveTpm() { return mRawEffectiveTpm; }
    public boolean tpmCalculated() { return mTpmCalculated; }

    public void setCalculatedTpms(double rawEffective, double effective, double expected)
    {
        mEffectiveTpm = effective;
        mRawEffectiveTpm = rawEffective;
        mExpectedTpm = expected;
        mTpmCalculated = true;
    }

    public String toString()
    {
        return format("%s neIds(%d) expected(%4.3e) effective(%4.3e raw=%4.3e)",
                Peptide, mNeIds.size(), mExpectedTpm, mEffectiveTpm, mRawEffectiveTpm);
    }
}
