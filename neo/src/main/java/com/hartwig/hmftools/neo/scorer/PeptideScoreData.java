package com.hartwig.hmftools.neo.scorer;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.PeptideData;
import com.hartwig.hmftools.neo.bind.BindData;

public class PeptideScoreData extends PeptideData
{
    private double mTpmUp;

    private final Map<Integer,Double> mNeoTpmDowns;
    private double mTpmDownTotal;

    private double mEffectiveTpm;
    private double mRawEffectiveTpm;

    private final List<BindData> mAlleleScoreData;

    public PeptideScoreData(final PeptideData peptideData, final int neId, final double tpmUp, final double tpmDown)
    {
        super(peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank);

        mTpmUp = tpmUp;
        mNeoTpmDowns = Maps.newHashMap();

        addNeoepitopeTpm(neId, tpmDown);

        mEffectiveTpm = 0;
        mRawEffectiveTpm = 0;

        mAlleleScoreData = Lists.newArrayList();
    }

    public void addAllele(final String allele)
    {
        BindData bindData = new BindData(allele, Peptide, "", UpFlank, DownFlank);
        bindData.setTPM(effectiveTpm());
        mAlleleScoreData.add(bindData);
    }

    public List<BindData> alleleScoreData() { return mAlleleScoreData; }

    public void addNeoepitopeTpm(final int neId, final double tpmDown)
    {
        mNeoTpmDowns.put(neId, tpmDown);
        mTpmDownTotal += tpmDown;
    }

    public boolean hasNeo(int neId) { return mNeoTpmDowns.containsKey(neId); }

    public double tpmUp() { return mTpmUp; }
    public double tpmDownTotal() { return mTpmDownTotal; }

    public double tpmUpAllocation(final int neId)
    {
        if(mNeoTpmDowns.size() == 1)
            return mTpmUp;

        Double neoTpmDown = mNeoTpmDowns.get(neId);
        if(neoTpmDown == null)
            return 0;

        return mTpmUp * neoTpmDown / mTpmDownTotal;
    }

    public double effectiveTpm() { return mEffectiveTpm; }
    public double rawEffectiveTpm() { return mRawEffectiveTpm; }

    public void setEffectiveTpms(double raw, double calc)
    {
        mEffectiveTpm = calc;
        mRawEffectiveTpm = raw;
    }

    public String toString()
    {
        return format("%s tpmUp(%4.3e) downTotal(%d=%4.3e) effective(%4.3e raw=%4.3e)",
                Peptide, mTpmUp, mNeoTpmDowns.size(), mTpmDownTotal, mEffectiveTpm, mRawEffectiveTpm);
    }
}
