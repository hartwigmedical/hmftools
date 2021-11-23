package com.hartwig.hmftools.isofox.unmapped;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.startEndStr;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class UnmappedRead
{
    public final String ReadId;
    public final ChrBaseRegion ReadRegion;
    public final int ScLength;
    public final int ScSide;
    public final double AvgBaseQual;
    public final String GeneName;
    public final String GeneId;
    public final String TransName;
    public final int ExonRank;
    public final int ExonBoundary;
    public final int ExonDistance;
    public final String SpliceType;
    public final String ScBases;
    public final String MateCoords;
    public final int CohortFrequency;
    public boolean MatchesSupplementary;

    private final String mPosKey;

    public static final String SPLICE_TYPE_ACCEPTOR = "acceptor";
    public static final String SPLICE_TYPE_DONOR = "donor";
    public static final String UMR_NO_MATE = "NO_MATE";

    public UnmappedRead(
            final String readId, final ChrBaseRegion readRegion, final int scLength, final int scSide,
            final double avgBaseQual, final String geneId, final String geneName, final String transName, final int exonRank,
            final int exonBoundary, final int exonDistance, final String spliceType, final String scBases, final String mateCoords,
            final int cohortFrequency, final boolean matchesSupplementary)
    {
        ReadId = readId;
        ReadRegion = readRegion;
        ScLength = scLength;
        ScSide = scSide;
        AvgBaseQual = avgBaseQual;
        GeneName = geneName;
        GeneId = geneId;
        TransName = transName;
        ExonRank = exonRank;
        ExonBoundary = exonBoundary;
        ExonDistance = exonDistance;
        SpliceType = spliceType;
        ScBases = scBases;
        MateCoords = mateCoords;
        CohortFrequency = cohortFrequency;
        MatchesSupplementary = matchesSupplementary;

        mPosKey = positionKey(ScSide, ExonBoundary);
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(DELIMITER);
        sj.add("ReadId");
        sj.add("Chromosome");
        sj.add("ReadStart");
        sj.add("ReadEnd");
        sj.add("SoftClipLength");
        sj.add("SoftClipSide");
        sj.add("SpliceType");
        sj.add("AvgBaseQual");
        sj.add("GeneId");
        sj.add("GeneName");
        sj.add("TransName");
        sj.add("ExonRank");
        sj.add("ExonBoundary");
        sj.add("ExonDistance");
        sj.add("SoftClipBases");
        sj.add("MateCoords");
        sj.add("CohortFreq");
        sj.add("MatchesSupp");

        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIMITER);
        sj.add(ReadId);
        sj.add(ReadRegion.Chromosome);
        sj.add(String.valueOf(ReadRegion.start()));
        sj.add(String.valueOf(ReadRegion.end()));
        sj.add(String.valueOf(ScLength));
        sj.add(String.valueOf(ScSide));
        sj.add(SpliceType);
        sj.add(String.format("%.1f", AvgBaseQual));
        sj.add(GeneId);
        sj.add(GeneName);
        sj.add(TransName);
        sj.add(String.valueOf(ExonRank));
        sj.add(String.valueOf(ExonBoundary));
        sj.add(String.valueOf(ExonDistance));
        sj.add(ScBases);
        sj.add(MateCoords);
        sj.add(String.valueOf(CohortFrequency));
        sj.add(String.valueOf(MatchesSupplementary));

        return sj.toString();
    }

    public String positionKey()
    {
        return mPosKey;
    }

    public static String positionKey(final int scSide, final int exonBoundary)
    {
        return String.format("%d_%d", exonBoundary, scSide);
    }

    public boolean matches(final UnmappedRead other)
    {
        return ReadRegion.Chromosome.equals(other.ReadRegion.Chromosome) && ExonBoundary == other.ExonBoundary && ScSide == other.ScSide;
    }

    public String toString()
    {
        return String.format("read(%s: %d -> %d) softClip(%s %s exon=%d dist=%d)",
                ReadRegion.Chromosome, ReadRegion.start(), ReadRegion.end(), startEndStr(ScSide), SpliceType, ExonBoundary, ExonDistance);
    }

}
