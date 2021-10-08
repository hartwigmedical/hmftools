package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.PaveConstants.DELIM;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.variant.Variant;

public class CodingContext
{
    public TranscriptRegionType RegionType; // favours more impactful type
    public TranscriptCodingType CodingType;
    public int ExonRank;

    public int CodingBase; // indexed at 1, coding base from start of transcript's coding region
    public int[] CodingPositionRange;
    public boolean SpansSpiceJunction;
    public boolean IsFrameShift;

    public int NonCodingBaseDistance;
    public int UpstreamPhase;
    public int BasesToLastExonJunction;

    public CodingContext()
    {
        RegionType = TranscriptRegionType.UNKNOWN;
        CodingType = TranscriptCodingType.UNKNOWN;
        ExonRank = 0;
        SpansSpiceJunction = false;
        IsFrameShift = false;
        CodingPositionRange = new int[] {0, 0};
        NonCodingBaseDistance = 0;
        UpstreamPhase = PHASE_NONE;
        BasesToLastExonJunction = -1;
    }

    public boolean isCoding() { return CodingType == CODING && RegionType == EXONIC; }

    public String codingRef(final VariantData variant)
    {
        if(variant.isInsert())
            return variant.Ref;

        return trimToCoding(variant.Ref, variant.Position, variant.isDeletion() ? variant.EndPosition - 1 : variant.EndPosition);
    }

    public String codingAlt(final VariantData variant)
    {
        if(variant.isInsert())
            return variant.Alt;

        if(variant.isDeletion())
        {
            // should it be trimmed to empty? yes if outside the coding region
            if(positionWithin(variant.Position, CodingPositionRange[SE_START], CodingPositionRange[SE_END]))
                return variant.Alt;
            else
                return "";
        }

        return trimToCoding(variant.Alt, variant.Position, variant.EndPosition);
    }

    private String trimToCoding(final String bases, int posStart, int posEnd)
    {
        if(CodingPositionRange[SE_START] > posStart)
            return bases.substring(CodingPositionRange[SE_START] - posStart);

        if(CodingPositionRange[SE_END] > 0 && CodingPositionRange[SE_END] < posEnd)
        {
            int diff = posEnd - CodingPositionRange[SE_END];
            return bases.substring(0, bases.length() - diff);
        }

        return bases;
    }

    public String hgvsStr() { return "tbc"; }

    public static String csvHeader()
    {
        return "HgvsCoding,RegionType,CodingType,ExonRank,CodingBase,CodingPosRange,UpstreamPhase,SpansSplice,NonCodingBaseDist";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIM);

        sj.add(hgvsStr());
        sj.add(RegionType.toString());
        sj.add(CodingType.toString());
        sj.add(String.valueOf(ExonRank));
        sj.add(String.valueOf(CodingBase));
        sj.add(CodingPositionRange[SE_START] == CodingPositionRange[SE_END] ? String.valueOf(CodingPositionRange[SE_START])
                : String.format("%d-%d", CodingPositionRange[SE_START], CodingPositionRange[SE_END]));
        sj.add(String.valueOf(UpstreamPhase));
        sj.add(String.valueOf(SpansSpiceJunction));
        sj.add(String.valueOf(NonCodingBaseDistance));

        return sj.toString();
    }
}
