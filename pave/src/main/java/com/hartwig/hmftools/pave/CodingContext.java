package com.hartwig.hmftools.pave;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcExonicCodingPhase;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.codingBaseLength;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.gene.CodingBaseData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;

// consider moving to hmf-common
public class CodingContext
{
    public TranscriptRegionType[] RegionTypes;
    public TranscriptCodingType[] CodingTypes;
    public int[] ExonRank;
    public int[] CodingBase; // index of coding base from start of transcript
    public int[] NonCodingBase; // bases before or after the nearest coding base (ie 5P, 3P or intronic)
    public int CodingPhase; // for upstream position
    public int BasesToLastExonJunction;

    public CodingContext()
    {
        RegionTypes = new TranscriptRegionType[] { TranscriptRegionType.UNKNOWN, TranscriptRegionType.UNKNOWN };
        CodingTypes = new TranscriptCodingType[] { TranscriptCodingType.UNKNOWN, TranscriptCodingType.UNKNOWN };
        ExonRank = new int[] {0, 0};
        CodingBase = new int[] {0, 0};
        NonCodingBase = new int[] {0, 0};
        CodingPhase = PHASE_NONE;
        BasesToLastExonJunction = -1;
    }

    public boolean isCoding()
    {
        return (CodingTypes[SE_START] == CODING && RegionTypes[SE_START] == EXONIC)
            || (CodingTypes[SE_END] == CODING && RegionTypes[SE_END] == EXONIC);
    }

    public static CodingContext determineContext(final VariantData variant, final TranscriptData transData)
    {
        int posStart = variant.Position;
        int posEnd = variant.EndPosition;

        if(posEnd < transData.TransStart || posStart > transData.TransEnd)
        {
            return determineUpstreamContext(posStart, posEnd, transData);
        }

        if(transData.CodingStart == null)
        {
            return determineNonCodingContext(posStart, posEnd, transData);
        }

        if(posStart < transData.CodingStart || posEnd > transData.CodingEnd)
        {
            return determinePreOrPostCodingContext(posStart, posEnd, transData);
        }

        // now handle scenario where variant is within the coding region
        // find the coding base index and if the variant is intronic, then the bases into the intron (closest to the exon)

        CodingContext cc = new CodingContext();

        int codingStart = transData.CodingStart;
        int codingEnd = transData.CodingEnd;
        int upstreamStartPos = transData.posStrand() ? variant.Position : variant.EndPosition;

        int preExonCodingBases = 0; // accumulated since start of coding region prior to exon closest to variant

        // method:
        // find intron or exon which contains the start position
        // if intronic then calculate distance to next exon
        // calculate coding bases prior to this location and add any addition if within an exon
        // set region types

        cc.CodingTypes[SE_START] = cc.CodingTypes[SE_END] = CODING;

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
                    int position = se == SE_START ? posStart : posEnd;

                    // accumulate coding base positions
                    if(positionWithin(position, exon.Start, exon.End))
                    {
                        cc.RegionTypes[se] = EXONIC;
                        cc.ExonRank[se] = exon.Rank;

                        // add in any extra coding bases
                        cc.CodingBase[se] = preExonCodingBases + position - max(codingStart, exon.Start) + 1;

                        cc.CodingPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, upstreamStartPos);
                    }
                    else if(nextExon != null && position > exon.End && position < nextExon.Start)
                    {
                        cc.CodingPhase = transData.Strand == POS_STRAND ? exon.PhaseEnd : nextExon.PhaseEnd;

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

                    continue;
                }

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    int position = se == SE_START ? posStart : posEnd;

                    // accumulate coding base positions
                    if(positionWithin(position, exon.Start, exon.End))
                    {
                        cc.RegionTypes[se] = EXONIC;
                        cc.ExonRank[se] = exon.Rank;

                        // add in any extra coding bases
                        cc.CodingBase[se] = preExonCodingBases + min(codingEnd, exon.End) - position + 1;

                        cc.CodingPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, upstreamStartPos);
                    }
                    else if(nextExon != null && position > nextExon.End && position < exon.Start)
                    {
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

    public static CodingContext determinePreOrPostCodingContext(final int posStart, final int posEnd, final TranscriptData transData)
    {
        CodingContext cc = new CodingContext();

        // record exonic bases from the position to the start of coding for 5'UTR, or end of coding to the position for 3'UTR
        if((transData.posStrand() && posStart > transData.CodingEnd) || (!transData.posStrand() && posEnd < transData.CodingStart))
        {
            // post-coding
            int codingBases = codingBaseLength(transData);
            cc.CodingBase[SE_START] = cc.CodingBase[SE_END] = codingBases;
            cc.CodingTypes[SE_START] = cc.CodingTypes[SE_END] = UTR_3P;
        }
        else
        {
            cc.CodingTypes[SE_START] = cc.CodingTypes[SE_END] = UTR_5P;
        }

        if(posStart < transData.CodingStart)
        {
            int codingStart = transData.CodingStart;

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                if(nextExon != null && posStart >= nextExon.Start)
                    continue;

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    int position = se == SE_START ? posStart : posEnd;

                    if(positionWithin(position, exon.Start, exon.End))
                    {
                        cc.ExonRank[se] = exon.Rank;
                        cc.RegionTypes[se] = EXONIC;

                        if(positionWithin(codingStart, exon.Start, exon.End))
                            cc.NonCodingBase[se] += codingStart - position;
                        else
                            cc.NonCodingBase[se] += exon.End - position;
                    }
                    else if(nextExon != null && position > exon.End && position < nextExon.Start)
                    {
                        cc.ExonRank[se] = transData.posStrand() ? exon.Rank : nextExon.Rank; // could be set to closest but not used for anything
                        cc.RegionTypes[se] = INTRONIC;
                    }
                    else if(exon.Start > position)
                    {
                        // past where the position is so now just about tracking the exonic bases before coding starts
                        if(positionWithin(codingStart, exon.Start, exon.End))
                            cc.NonCodingBase[se] += codingStart - exon.Start;
                        else
                            cc.NonCodingBase[se] += exon.baseLength();
                    }
                }

                if(positionWithin(codingStart, exon.Start, exon.End))
                    break;
            }
        }
        else
        {
            // similar for position after coding end
            int codingEnd = transData.CodingEnd;

            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i >= 1 ? transData.exons().get(i - 1) : null;

                if(nextExon != null && posEnd <= nextExon.End)
                    continue;

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    int position = se == SE_START ? posStart : posEnd;

                    if(positionWithin(position, exon.Start, exon.End))
                    {
                        cc.ExonRank[se] = exon.Rank;
                        cc.RegionTypes[se] = EXONIC;

                        if(positionWithin(codingEnd, exon.Start, exon.End))
                            cc.NonCodingBase[se] += position - codingEnd;
                        else
                            cc.NonCodingBase[se] += position - exon.Start;
                    }
                    else if(nextExon != null && position > nextExon.End && position < exon.Start)
                    {
                        cc.ExonRank[se] = transData.posStrand() ? nextExon.Rank : exon.Rank;
                        cc.RegionTypes[se] = INTRONIC;
                    }
                    else if(exon.End < position)
                    {
                        if(positionWithin(codingEnd, exon.Start, exon.End))
                            cc.NonCodingBase[se] += exon.End - codingEnd;
                        else
                            cc.NonCodingBase[se] += exon.baseLength();
                    }
                }

                if(positionWithin(codingEnd, exon.Start, exon.End))
                    break;
            }
        }

        // convention is to have negative values if 5'UTR
        if(cc.CodingTypes[SE_START] == UTR_5P)
        {
            cc.NonCodingBase[SE_START] *= -1;
            cc.NonCodingBase[SE_END] *= -1;
        }

        return cc;
    }

    public static CodingContext determineNonCodingContext(final int posStart, final int posEnd, final TranscriptData transData)
    {
        CodingContext cc = new CodingContext();

        // can set exon ranks and bases from start of transcript but that's it
        cc.CodingTypes[SE_START] = NON_CODING;
        cc.CodingTypes[SE_END] = NON_CODING;

        // how to define bases for HGVC coding if at all?
        if(transData.posStrand())
        {
            cc.NonCodingBase[SE_START] = posStart - transData.TransStart;
            cc.NonCodingBase[SE_END] = posEnd - transData.TransStart;

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);

                if(posStart > exon.End)
                    continue;

                if(posStart < exon.Start)
                {
                    cc.ExonRank[SE_START] = exon.Rank - 1;
                    cc.RegionTypes[SE_START] = INTRONIC;
                }
                else if(positionWithin(posStart, exon.Start, exon.End))
                {
                    cc.ExonRank[SE_START] = exon.Rank;
                    cc.RegionTypes[SE_START] = EXONIC;
                }

                if(posEnd < exon.Start)
                {
                    cc.ExonRank[SE_END] = exon.Rank - 1;
                    cc.RegionTypes[SE_END] = INTRONIC;
                }
                else if(positionWithin(posEnd, exon.Start, exon.End))
                {
                    cc.ExonRank[SE_END] = exon.Rank;
                    cc.RegionTypes[SE_END] = EXONIC;
                }

                break;
            }
        }
        else
        {
            cc.NonCodingBase[SE_START] = transData.TransEnd - posStart;
            cc.NonCodingBase[SE_END] = transData.TransEnd - posEnd;

            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                final ExonData exon = transData.exons().get(i);

                if(posEnd < exon.Start)
                    continue;

                if(posEnd > exon.End)
                {
                    cc.ExonRank[SE_END] = exon.Rank - 1;
                    cc.RegionTypes[SE_END] = INTRONIC;
                }
                else if(positionWithin(posEnd, exon.Start, exon.End))
                {
                    cc.ExonRank[SE_END] = exon.Rank;
                    cc.RegionTypes[SE_END] = EXONIC;
                }

                if(posStart > exon.End)
                {
                    cc.ExonRank[SE_START] = exon.Rank - 1;
                    cc.RegionTypes[SE_START] = INTRONIC;
                }
                else if(positionWithin(posStart, exon.Start, exon.End))
                {
                    cc.ExonRank[SE_START] = exon.Rank;
                    cc.RegionTypes[SE_START] = EXONIC;
                }

                break;
            }
        }

        return cc;
    }

    public static CodingContext determineUpstreamContext(final int posStart, final int posEnd, final TranscriptData transData)
    {
        CodingContext cc = new CodingContext();

        cc.RegionTypes[SE_START] = UPSTREAM;
        cc.RegionTypes[SE_END] = UPSTREAM;

        // upstream (or downstream but not handled)
        if(transData.posStrand())
        {
            // will be negative
            cc.NonCodingBase[SE_START] = posStart - transData.TransStart;
            cc.NonCodingBase[SE_END] = posEnd - transData.TransStart;
        }
        else
        {
            cc.NonCodingBase[SE_START] = transData.TransEnd - posStart;
            cc.NonCodingBase[SE_END] = transData.TransEnd - posEnd;
        }

        return cc;
    }
}
