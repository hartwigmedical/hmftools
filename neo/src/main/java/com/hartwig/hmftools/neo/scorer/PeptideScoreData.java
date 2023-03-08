package com.hartwig.hmftools.neo.scorer;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.PeptideData;
import com.hartwig.hmftools.neo.bind.BindData;

public class PeptideScoreData extends PeptideData
{
    private double mTpmUp;

    private final Map<Integer,Double> mNeoTpmDowns;
    private double mTpmDownTotal;

    private double mExpectedTpm;
    private double mEffectiveTpm;
    private double mRawEffectiveTpm;

    private final List<BindData> mAlleleScoreData;

    private boolean mWritten;

    public PeptideScoreData(final PeptideData peptideData, final int neId, final double tpmUp, final double tpmDown)
    {
        super(peptideData.Peptide, peptideData.UpFlank, peptideData.DownFlank);

        mTpmUp = tpmUp;
        mNeoTpmDowns = Maps.newHashMap();

        addNeoepitopeTpm(neId, tpmDown);

        mExpectedTpm = 0;
        mEffectiveTpm = 0;
        mRawEffectiveTpm = 0;

        mAlleleScoreData = Lists.newArrayList();
        mWritten = false;
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

    public boolean written() { return mWritten; }
    public void setWritten() { mWritten = true; }

    public boolean hasNeo(int neId) { return mNeoTpmDowns.containsKey(neId); }

    public String neoIdsStr()
    {
        if(mNeoTpmDowns.size() == 1)
            return String.valueOf(mNeoTpmDowns.keySet().iterator().next());

        StringJoiner sj = new StringJoiner(";");
        mNeoTpmDowns.keySet().forEach(x -> sj.add(String.valueOf(x)));
        return sj.toString();
    }

    public double tpmUp() { return mTpmUp; }
    public double tpmDownTotal() { return mTpmDownTotal; }

    public double expectedTpm() { return mExpectedTpm; }
    public double effectiveTpm() { return mEffectiveTpm; }
    public double rawEffectiveTpm() { return mRawEffectiveTpm; }

    public void setCalculatedTpms(double rawEffective, double effective, double expected)
    {
        mEffectiveTpm = effective;
        mRawEffectiveTpm = rawEffective;
        mExpectedTpm = expected;
    }

    public String toString()
    {
        return format("%s tpmUp(%4.3e) downTotal(%d=%4.3e) expected(%4.3e) effective(%4.3e raw=%4.3e)",
                Peptide, mTpmUp, mNeoTpmDowns.size(), mTpmDownTotal, mExpectedTpm, mEffectiveTpm, mRawEffectiveTpm);
    }
}
