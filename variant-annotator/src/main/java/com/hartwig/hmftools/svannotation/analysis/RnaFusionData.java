package com.hartwig.hmftools.svannotation.analysis;

public class RnaFusionData
{
    public final String Name;
    public final String GeneUp;
    public final String GeneDown;
    public final long PositionUp;
    public final long PositionDown;
    public final byte OrientUp;
    public final byte OrientDown;

    private int mExonMinRankUp;
    private int mExonMaxRankUp;
    private int mExonMinRankDown;
    private int mExonMaxRankDown;

    public RnaFusionData(final String name, final String geneUp, final String geneDown,
            long posUp, long posDown, byte orientUp, byte orientDown)
    {
        Name = name;
        GeneUp = geneUp;
        GeneDown = geneDown;
        PositionUp = posUp;
        PositionDown = posDown;
        OrientUp = orientUp;
        OrientDown = orientDown;

        mExonMinRankUp = 0;
        mExonMaxRankUp = 0;
        mExonMinRankDown = 0;
        mExonMaxRankDown = 0;

    }

    public void setExonUpRank(int min, int max)
    {
        mExonMaxRankUp = max;
        mExonMinRankUp = min;
    }

    public void setExonDownRank(int min, int max)
    {
        mExonMaxRankDown = max;
        mExonMinRankDown = min;
    }

    public int exonMinRankUp() { return mExonMinRankUp; }
    public int exonMaxRankUp() {return mExonMaxRankUp; }
    public int exonMinRankDown() { return mExonMinRankDown; }
    public int exonMaxRankDown() { return mExonMaxRankDown; }

}
