package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.CodingUtils.determineCodingContext;
import static com.hartwig.hmftools.pave.NonCodingContext.determineNonCodingContext;
import static com.hartwig.hmftools.pave.NonCodingContext.determinePreOrPostCodingContext;
import static com.hartwig.hmftools.pave.NonCodingContext.determineUpstreamContext;
import static com.hartwig.hmftools.pave.PaveConstants.DELIM;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;

public class CodingContext
{
    public TranscriptRegionType RegionType; // favours more impactful type
    public TranscriptCodingType CodingType;
    public int ExonRank;

    public int CodingBase; // indexed at 1, coding base from start of transcript's coding region
    public int[] CodingPositionRange;
    public boolean SpansSpiceJunction;
    public int DeletedCodingBases;

    public int NonCodingBaseDistance;
    public int UpstreamPhase;
    public int BasesToLastExonJunction;

    public CodingContext()
    {
        RegionType = TranscriptRegionType.UNKNOWN;
        CodingType = TranscriptCodingType.UNKNOWN;
        ExonRank = 0;
        SpansSpiceJunction = false;
        CodingPositionRange = new int[] {0, 0};
        DeletedCodingBases = 0;
        NonCodingBaseDistance = 0;
        UpstreamPhase = PHASE_NONE;
        BasesToLastExonJunction = -1;
    }

    public boolean isCoding() { return CodingType == CODING && RegionType == EXONIC; }

    public String hgvsStr() { return "tbc"; }

    public static CodingContext determineContext(final VariantData variant, final TranscriptData transData)
    {
        int posStart = variant.Position;
        int posEnd = variant.EndPosition;

        if(variant.altBasesBelow(transData.TransStart) || variant.altBasesAbove(transData.TransEnd))
        {
            return determineUpstreamContext(posStart, posEnd, transData);
        }

        if(transData.nonCoding())
        {
            return determineNonCodingContext(posStart, posEnd, transData);
        }

        if(variant.altBasesBelow(transData.CodingStart) || variant.altBasesAbove(transData.CodingEnd))
        {
            return determinePreOrPostCodingContext(variant, transData);
        }

        return determineCodingContext(variant, transData);
    }

    /*
    private static CodingContext determineCodingContext(final VariantData variant, final TranscriptData transData)
    {
        // now handle scenario where variant is within the coding region
        // find the coding base index and if the variant is intronic, then the bases into the intron (closest to the exon)

        CodingContext cc = new CodingContext();

        int posStart = variant.Position;
        int posEnd = variant.EndPosition;
        int codingStart = transData.CodingStart;
        int codingEnd = transData.CodingEnd;
        int upstreamStartPos = transData.posStrand() ? variant.Position : variant.EndPosition;

        int preExonCodingBases = 0; // accumulated since start of coding region prior to exon closest to variant

        // method:
        // find intron or exon which contains the start position
        // if intronic then calculate distance to next exon
        // calculate coding bases prior to this location and add any addition if within an exon
        // set region types

        cc.CodingType[SE_START] = cc.CodingType[SE_END] = CODING;

        if(transData.posStrand())
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                if(posStart > exon.End)
                {
                    // keep track of coding bases if this exon overlaps the coding region
                    if(positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        preExonCodingBases += min(codingEnd, exon.End) - max(codingStart, exon.Start);

                    if(nextExon != null && posStart >= nextExon.Start)
                        continue;
                }

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    if(se == SE_END && posStart == posEnd)
                    {
                        cc.copyStartValues();
                        continue;
                    }

                    int position = se == SE_START ? posStart : posEnd;

                    // accumulate coding base positions
                    if(positionWithin(position, exon.Start, exon.End))
                    {
                        cc.RegionType[se] = EXONIC;
                        cc.ExonRank[se] = exon.Rank;

                        // add in any extra coding bases
                        cc.CodingBase[se] = preExonCodingBases + position - max(codingStart, exon.Start) + 1;

                        cc.CodingPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, upstreamStartPos);
                    }
                    else if(nextExon != null && position > exon.End && position < nextExon.Start)
                    {
                        cc.CodingPhase = transData.Strand == POS_STRAND ? exon.PhaseEnd : nextExon.PhaseEnd;
                        cc.RegionType[se] = INTRONIC;

                        int distanceToPrev = position - exon.End;
                        int distanceToNext = nextExon.Start - position;

                        if(distanceToPrev < distanceToNext)
                        {
                            cc.ExonRank[se] = exon.Rank;
                            cc.CodingBase[se] = preExonCodingBases;
                            cc.NonCodingBase[se] = position - exon.End;
                        }
                        else
                        {
                            cc.ExonRank[se] = nextExon.Rank;
                            cc.CodingBase[se] = preExonCodingBases + 1; // first base of next exon
                            cc.NonCodingBase[se] = position - nextExon.Start; // -ve for back into previous intron
                        }
                    }
                }

                break;
            }
        }
        else
        {
            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i >= 1 ? transData.exons().get(i - 1) : null;

                if(posEnd < exon.Start)
                {
                    if(positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        preExonCodingBases += min(codingEnd, exon.End) - max(codingStart, exon.Start);

                    if(nextExon != null && posEnd <= nextExon.End)
                        continue;
                }

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    if(se == SE_END && posStart == posEnd)
                    {
                        cc.copyStartValues();
                        continue;
                    }

                    int position = se == SE_START ? posStart : posEnd;

                    // accumulate coding base positions
                    if(positionWithin(position, exon.Start, exon.End))
                    {
                        cc.RegionType[se] = EXONIC;
                        cc.ExonRank[se] = exon.Rank;

                        // add in any extra coding bases
                        cc.CodingBase[se] = preExonCodingBases + min(codingEnd, exon.End) - position + 1;

                        cc.CodingPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, upstreamStartPos);
                    }
                    else if(nextExon != null && position > nextExon.End && position < exon.Start)
                    {
                        cc.RegionType[se] = INTRONIC;
                        cc.CodingPhase = transData.Strand == POS_STRAND ? exon.PhaseEnd : nextExon.PhaseEnd;

                        int distanceToPrev = exon.Start - position;
                        int distanceToNext = position - nextExon.End;

                        if(distanceToPrev < distanceToNext)
                        {
                            cc.ExonRank[se] = exon.Rank;
                            cc.CodingBase[se] = preExonCodingBases;
                            cc.NonCodingBase[se] = exon.Start - position;
                        }
                        else
                        {
                            cc.ExonRank[se] = nextExon.Rank;
                            cc.CodingBase[se] = preExonCodingBases + 1; // first base of next exon
                            cc.NonCodingBase[se] = nextExon.End - position; // negative
                        }
                    }
                }

                break;
            }
        }

        return cc;
    }
    */

    public static String csvHeader()
    {
        return "HgvsCoding,RegionType,CodingType,ExonRank,CodingBase,CodingPosRange,SpansSplice,NonCodingBaseDist";
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
        sj.add(String.valueOf(SpansSpiceJunction));
        sj.add(String.valueOf(NonCodingBaseDistance));

        return sj.toString();
    }
}
