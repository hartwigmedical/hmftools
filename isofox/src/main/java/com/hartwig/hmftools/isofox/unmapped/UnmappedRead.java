package com.hartwig.hmftools.isofox.unmapped;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.startEndStr;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

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
    public final byte[] ScBasesQuals;
    public final String MateCoords;
    public boolean MatchesChimeric;

    private final String mPosKey;

    public static final String SPLICE_TYPE_ACCEPTOR = "acceptor";
    public static final String SPLICE_TYPE_DONOR = "donor";
    public static final String UMR_NO_MATE = "NO_MATE";

    public static final int MAX_BASE_QUAL = 36;
    public static final char MAX_BASE_QUAL_ID = 'm';

    public UnmappedRead(
            final String readId, final ChrBaseRegion readRegion, final int scLength, final int scSide,
            final double avgBaseQual, final String geneId, final String geneName, final String transName, final int exonRank,
            final int exonBoundary, final int exonDistance, final String spliceType, final String scBases, final byte[] scBasesQuals,
            final String mateCoords, final boolean matchesChimeric)
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
        ScBasesQuals = scBasesQuals;
        MateCoords = mateCoords;
        MatchesChimeric = matchesChimeric;

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
        sj.add("MateCoords");
        sj.add("MatchesChimeric");
        sj.add("SoftClipBases");
        sj.add("SoftClipBaseQuals");

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
        sj.add(MateCoords);
        sj.add(String.valueOf(MatchesChimeric));
        sj.add(ScBases);
        sj.add(baseQualsToString(ScBasesQuals));

        return sj.toString();
    }

    public static byte[] baseQualsFromString(final String qualString, int baseCount)
    {
        byte[] quals = new byte[baseCount];

        int index = 0;

        for(int i = 0; i < baseCount; ++i)
        {
            char qual = qualString.charAt(index);

            if(qual == MAX_BASE_QUAL_ID)
            {
                quals[i] = MAX_BASE_QUAL;
                ++index;
            }
            else
            {
                quals[i] = Byte.parseByte(qualString.substring(index, index + 2));
                index += 2;
            }
        }

        return quals;
    }

    public static String baseQualsToString(final byte[] quals)
    {
        if(quals == null || quals.length == 0)
            return "";

        StringJoiner sj = new StringJoiner("");

        for(int i = 0; i < quals.length; ++i)
        {
            if(quals[i] == MAX_BASE_QUAL)
                sj.add(String.valueOf(MAX_BASE_QUAL_ID));
            else
                sj.add(String.format("%02d", quals[i]));
        }

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
