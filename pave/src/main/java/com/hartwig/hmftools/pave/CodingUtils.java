package com.hartwig.hmftools.pave;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcExonicCodingPhase;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.codingBaseLength;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public final class CodingUtils
{
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

    public static CodingContext determineCodingContext(final VariantData variant, final TranscriptData transData)
    {
        // initially set the region for this variant and the coding base range and phase if any
        // calculate coding bases prior to the upstream position and add any addition if within an exon

        CodingContext cc = new CodingContext();

        cc.CodingType = CODING;

        int codingStart = transData.CodingStart;
        int codingEnd = transData.CodingEnd;

        boolean posStrand = transData.posStrand();
        int preExonCodingBases = 0; // accumulated since start of coding region prior to exon closest to variant

        int upstreamStartPos = posStrand ? variant.Position : variant.EndPosition;

        // INDELs don't affect the last base of an exon so need diff criteria to skip it
        int skipExonPos = (!variant.isIndel() == posStrand) ? variant.Position : variant.EndPosition;

        if(transData.posStrand())
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                if(skipExonPos > exon.End)
                {
                    // keep track of coding bases if this exon overlaps the coding region but is before the variant
                    if(positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        preExonCodingBases += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

                    if(nextExon != null && upstreamStartPos >= nextExon.Start)
                        continue;
                }

                // check for fully intronic, fully exonic or spanning
                boolean withinExon = variant.altPositionsWithin(exon.Start, exon.End);
                if(withinExon || variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    cc.RegionType = EXONIC;
                    cc.ExonRank = exon.Rank;
                    cc.SpansSpiceJunction = !withinExon;

                    cc.UpstreamPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, max(upstreamStartPos, exon.Start));

                    cc.CodingPositionRange[SE_START] = max(max(exon.Start, codingStart), variant.Position);

                    cc.CodingPositionRange[SE_END] = !variant.isInsert() ?
                            min(min(exon.End, codingEnd), variant.EndPosition) : cc.CodingPositionRange[SE_START];

                    // add in any extra coding bases
                    cc.CodingBase = preExonCodingBases + cc.CodingPositionRange[SE_START] - max(codingStart, exon.Start) + 1;
                }
                else if(nextExon != null && variant.altPositionsOverlap(nextExon.Start, nextExon.End))
                {
                    continue;
                }
                else
                {
                    cc.UpstreamPhase = posStrand ? exon.PhaseEnd : nextExon.PhaseEnd;
                    cc.RegionType = INTRONIC;

                    int distanceToPrev = upstreamStartPos - exon.End;
                    int distanceToNext = nextExon.Start - upstreamStartPos;

                    if(distanceToPrev < distanceToNext)
                    {
                        cc.ExonRank = exon.Rank;
                        cc.CodingBase = preExonCodingBases;
                        cc.NonCodingBaseDistance = upstreamStartPos - exon.End;
                    }
                    else
                    {
                        cc.ExonRank = nextExon.Rank;
                        cc.CodingBase = preExonCodingBases + 1; // first base of next exon
                        cc.NonCodingBaseDistance = upstreamStartPos - nextExon.Start; // -ve for back into previous intron
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

                if(skipExonPos < exon.Start)
                {
                    if(positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        preExonCodingBases += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

                    if(nextExon != null && upstreamStartPos <= nextExon.End)
                        continue;
                }

                boolean withinExon = variant.altPositionsWithin(exon.Start, exon.End);
                if(withinExon || variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    cc.RegionType = EXONIC;
                    cc.ExonRank = exon.Rank;
                    cc.SpansSpiceJunction = !withinExon;

                    cc.UpstreamPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, min(upstreamStartPos, exon.End));

                    cc.CodingPositionRange[SE_START] = max(max(exon.Start, codingStart), variant.Position);
                    cc.CodingPositionRange[SE_END] = !variant.isInsert() ?
                            min(min(exon.End, codingEnd), variant.EndPosition) : cc.CodingPositionRange[SE_START];

                    // add in any extra coding bases
                    cc.CodingBase = preExonCodingBases + min(codingEnd, exon.End) - cc.CodingPositionRange[SE_END] + 1;
                }
                else if(nextExon != null && variant.altPositionsOverlap(nextExon.Start, nextExon.End))
                {
                    continue;
                }
                else
                {
                    cc.UpstreamPhase = posStrand ? exon.PhaseEnd : nextExon.PhaseEnd;
                    cc.RegionType = INTRONIC;

                    int distanceToPrev = exon.Start - upstreamStartPos;
                    int distanceToNext = upstreamStartPos - nextExon.End;

                    if(distanceToPrev < distanceToNext)
                    {
                        cc.ExonRank = exon.Rank;
                        cc.CodingBase = preExonCodingBases;
                        cc.NonCodingBaseDistance = exon.Start - upstreamStartPos;
                    }
                    else
                    {
                        cc.ExonRank = nextExon.Rank;
                        cc.CodingBase = preExonCodingBases + 1; // first base of next exon
                        cc.NonCodingBaseDistance = nextExon.End - upstreamStartPos; // negative
                    }
                }

                break;
            }
        }

        if(variant.isIndel() && cc.CodingType == CODING && cc.RegionType == EXONIC)
        {
            // now coding bases have been clarified, mark any frameshift INDELs
            int adjustedCodingBases;

            if(variant.isInsert())
            {
                adjustedCodingBases = variant.baseDiff();
            }
            else
            {
                String ref = cc.codingRef(variant);
                String alt = cc.codingAlt(variant);

                adjustedCodingBases = variant.isInsert() ? variant.baseDiff() : ref.length() - alt.length();
            }

            if((adjustedCodingBases % 3) != 0)
                cc.IsFrameShift = true;
        }

        return cc;
    }

    public static CodingContext determinePreOrPostCodingContext(final VariantData variant, final TranscriptData transData)
    {
        CodingContext cc = new CodingContext();

        // record exonic bases from the position to the start of coding for 5'UTR, or end of coding to the position for 3'UTR
        if((transData.posStrand() && variant.altBasesAbove(transData.CodingEnd))
                || (!transData.posStrand() && variant.altBasesBelow(transData.CodingStart)))
        {
            // post-coding
            cc.CodingBase = codingBaseLength(transData);
            cc.CodingType = UTR_3P;
        }
        else
        {
            cc.CodingType = UTR_5P;
        }

        if(variant.altBasesBelow(transData.CodingStart))
        {
            int codingStart = transData.CodingStart;
            int position = variant.EndPosition; // TO-DO - use alt methods where possible

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                if(nextExon != null && variant.altBasesAbove(nextExon.Start - 1))
                    continue;

                if(variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    cc.ExonRank = exon.Rank;
                    cc.RegionType = EXONIC;

                    if(positionWithin(codingStart, exon.Start, exon.End))
                        cc.NonCodingBaseDistance += codingStart - position;
                    else
                        cc.NonCodingBaseDistance += exon.End - position;
                }
                else if(nextExon != null && position > exon.End && position < nextExon.Start)
                {
                    cc.ExonRank = transData.posStrand() ? exon.Rank : nextExon.Rank; // could be set to closest but not used for anything
                    cc.RegionType = INTRONIC;
                }
                else if(exon.Start > position)
                {
                    // past where the position is so now just about tracking the exonic bases before coding starts
                    if(positionWithin(codingStart, exon.Start, exon.End))
                        cc.NonCodingBaseDistance += codingStart - exon.Start;
                    else
                        cc.NonCodingBaseDistance += exon.baseLength();
                }

                if(positionWithin(codingStart, exon.Start, exon.End))
                    break;
            }
        }
        else
        {
            // similar for position after coding end
            int codingEnd = transData.CodingEnd;
            int position = variant.Position; // TO-DO - use alt methods where possible

            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i >= 1 ? transData.exons().get(i - 1) : null;

                if(nextExon != null && variant.altBasesBelow(nextExon.End + 1))
                    continue;

                if(variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    cc.ExonRank = exon.Rank;
                    cc.RegionType = EXONIC;

                    if(positionWithin(codingEnd, exon.Start, exon.End))
                        cc.NonCodingBaseDistance += position - codingEnd;
                    else
                        cc.NonCodingBaseDistance += position - exon.Start;
                }
                else if(nextExon != null && position > nextExon.End && position < exon.Start)
                {
                    cc.ExonRank = transData.posStrand() ? nextExon.Rank : exon.Rank;
                    cc.RegionType = INTRONIC;
                }
                else if(exon.End < position)
                {
                    if(positionWithin(codingEnd, exon.Start, exon.End))
                        cc.NonCodingBaseDistance += exon.End - codingEnd;
                    else
                        cc.NonCodingBaseDistance += exon.baseLength();
                }

                if(positionWithin(codingEnd, exon.Start, exon.End))
                    break;
            }
        }

        // convention is to have negative values if 5'UTR
        if(cc.CodingType == UTR_5P)
        {
            cc.NonCodingBaseDistance *= -1;
        }

        return cc;
    }

    public static CodingContext determineNonCodingContext(final int posStart, final int posEnd, final TranscriptData transData)
    {
        CodingContext cc = new CodingContext();

        // can set exon ranks and bases from start of transcript but that's it
        cc.CodingType = NON_CODING;

        // how to define bases for HGVC coding if at all?
        if(transData.posStrand())
        {
            cc.NonCodingBaseDistance = posStart - transData.TransStart;

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);

                if(posStart > exon.End)
                    continue;

                if(posStart < exon.Start)
                {
                    cc.ExonRank = exon.Rank - 1;
                    cc.RegionType = INTRONIC;
                }
                else if(positionWithin(posStart, exon.Start, exon.End))
                {
                    cc.ExonRank = exon.Rank;
                    cc.RegionType = EXONIC;
                }

                break;
            }
        }
        else
        {
            cc.NonCodingBaseDistance = transData.TransEnd - posEnd;

            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                final ExonData exon = transData.exons().get(i);

                if(posEnd < exon.Start)
                    continue;

                if(posEnd > exon.End)
                {
                    cc.ExonRank = exon.Rank - 1;
                    cc.RegionType = INTRONIC;
                }
                else if(positionWithin(posEnd, exon.Start, exon.End))
                {
                    cc.ExonRank = exon.Rank;
                    cc.RegionType = EXONIC;
                }

                break;
            }
        }

        return cc;
    }

    public static CodingContext determineUpstreamContext(final int posStart, final int posEnd, final TranscriptData transData)
    {
        CodingContext cc = new CodingContext();

        cc.RegionType = UPSTREAM;

        // upstream (or downstream but not handled)
        if(transData.posStrand())
        {
            // will be negative
            cc.NonCodingBaseDistance = posStart - transData.TransStart;
        }
        else
        {
            cc.NonCodingBaseDistance = transData.TransEnd - posStart;
        }

        return cc;
    }

}
