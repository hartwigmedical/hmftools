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

public class CodingContext
{
    public TranscriptRegionType RegionType; // favours more impactful type if it spans
    public TranscriptCodingType CodingType;
    public int ExonRank;
    public byte Strand; // for convenience
    public int UpstreamPhase;

    // coding base from start of transcript's coding region, or exonic bases to coding if UTR, or exonic index if non-coding
    //
    public int CodingBase; // indexed at 1, for non-coding and UTR regions, this is exonic bases from the coding region

    // the alt bases for an MNV, otherwise the start and end ref bases (ie enclosing the del or insert)
    // likewise the exonic base range for non-coding and UTR regions
    public int[] CodingPositionRange;

    public int NearestExonDistance;
    public boolean SpansSpliceJunction;
    public boolean IsFrameShift;
    public boolean SpansCodingStart;
    public boolean SpansCodingEnd;
    public int DeletedCodingBases;

    public String SpliceDonorAcceptorBases;
    public boolean CodingEndsOnExonBoundary; // applicable for UTR if position is in the intron immediate before or after

    public String Hgvs;

    public CodingContext()
    {
        RegionType = TranscriptRegionType.UNKNOWN;
        CodingType = TranscriptCodingType.UNKNOWN;
        ExonRank = 0;
        Strand = 0;

        CodingBase = 0;
        CodingPositionRange = new int[] {0, 0};
        SpansSpliceJunction = false;
        IsFrameShift = false;
        SpliceDonorAcceptorBases = "";
        DeletedCodingBases = 0;

        NearestExonDistance = 0;
        UpstreamPhase = PHASE_NONE;
        CodingEndsOnExonBoundary = false;
        SpansCodingStart = false;
        SpansCodingEnd = false;
        Hgvs = "";
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
        int startDiff = CodingPositionRange[SE_START] > posStart ? CodingPositionRange[SE_START] - posStart : 0;
        int endDiff = CodingPositionRange[SE_END] > 0 && CodingPositionRange[SE_END] < posEnd ? posEnd - CodingPositionRange[SE_END] : 0;
        return bases.substring(startDiff, bases.length() - endDiff);
    }

    public static String csvHeader()
    {
        return "HgvsCoding,RegionType,CodingType,ExonRank,CodingBase,CodingPosRange,UpstreamPhase,SpansSplice,SpliceDonorAcceptorBases,NearestExonDistance";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIM);

        sj.add(Hgvs);
        sj.add(RegionType.toString());
        sj.add(CodingType.toString());
        sj.add(String.valueOf(ExonRank));
        sj.add(String.valueOf(CodingBase));
        sj.add(CodingPositionRange[SE_START] == CodingPositionRange[SE_END] ? String.valueOf(CodingPositionRange[SE_START])
                : String.format("%d-%d", CodingPositionRange[SE_START], CodingPositionRange[SE_END]));
        sj.add(String.valueOf(UpstreamPhase));
        sj.add(String.valueOf(SpansSpliceJunction));
        sj.add(SpliceDonorAcceptorBases);
        sj.add(String.valueOf(NearestExonDistance));

        return sj.toString();
    }
}
