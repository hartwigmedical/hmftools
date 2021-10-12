package com.hartwig.hmftools.pave;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
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

public final class CodingUtils
{
    public static void determineContext(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        cc.Strand = transData.Strand;

        setCodingContext(variant, transData, cc);

        // set pre-coding exonic bases separately since doesn't fit with way exons are iterated through
        if(cc.CodingType == UTR_5P)
            setPreCodingExonicDistance(variant, transData, cc);

        /*
        if(variant.altBasesBelow(transData.TransStart) || variant.altBasesAbove(transData.TransEnd))
        {
            setUpstreamContext(posStart, posEnd, transData, cc);
        }
        else if(transData.nonCoding())
        {
            return setNonCodingContext(posStart, posEnd, transData, cc);
        }
        else if(variant.altBasesBelow(transData.CodingStart) || variant.altBasesAbove(transData.CodingEnd))
        {
            return setPreOrPostCodingContext(variant, transData, cc);
        }
        else
        {
            setCodingContext(variant, transData, cc);
        }
        */
    }

    private static void setCodingContext(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        // set the following data
        // region type
        // coding type
        // distance to nearest exon (or upstream of transcript start)
        // coding base if within coding or number of exonic bases to coding start or end if 5' or 3' UTR
        //
        // other coding fields

        boolean inCoding = false;
        boolean isCodingTrans = !transData.nonCoding();
        Integer codingStart = transData.CodingStart;
        Integer codingEnd = transData.CodingEnd;

        if(variant.altBasesBelow(transData.TransStart) || variant.altBasesAbove(transData.TransEnd))
        {
            setUpstreamContext(variant, transData, cc);
            return;
        }

        boolean posStrand = transData.posStrand();

        if(transData.nonCoding())
        {
            cc.CodingType = NON_CODING;
        }
        else
        {
            boolean belowCodingStart = variant.altBasesBelow(codingStart);
            boolean aboveCodingEnd = variant.altBasesAbove(codingEnd);

            if(belowCodingStart)
            {
                cc.CodingType = posStrand ? UTR_5P : UTR_3P;
            }
            else if(aboveCodingEnd)
            {
                cc.CodingType = posStrand ? UTR_3P : UTR_5P;
            }
            else
            {
                inCoding = true;
                cc.CodingType = CODING;
            }
        }

        int prePosExonicBaseCount = 0; // will be used for non-coding transcripts only
        int postCodingExonicBases = 0; // exonic bases from the end of coding to the position for 3'UTR
        int codingBaseCount = 0; // accumulated since start of coding region prior to exon closest to variant

        int upstreamStartPos = posStrand ? variant.Position : variant.EndPosition;

        // INDELs don't affect the last base of an exon so need diff criteria to skip it
        int skipExonPos = (!variant.isIndel() == posStrand) ? variant.Position : variant.EndPosition;

        if(posStrand)
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                if(skipExonPos > exon.End)
                {
                    // keep track of coding bases if this exon overlaps the coding region but is before the variant
                    if(inCoding && positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        codingBaseCount += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

                    prePosExonicBaseCount += exon.baseLength();

                    if(isCodingTrans && exon.End >= codingEnd)
                        postCodingExonicBases += exon.End - max(codingEnd + 1, exon.Start) + 1;

                    if(nextExon != null && upstreamStartPos >= nextExon.Start)
                        continue;
                }

                // check for fully intronic, fully exonic or spanning
                boolean withinExon = variant.altPositionsWithin(exon.Start, exon.End);

                if(withinExon || variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    cc.RegionType = EXONIC;
                    cc.ExonRank = exon.Rank;
                    cc.SpansSpliceJunction = !withinExon;

                    prePosExonicBaseCount += upstreamStartPos - exon.Start + 1;

                    // bases from end of coding to the position
                    if(isCodingTrans && exon.End >= codingEnd)
                        postCodingExonicBases += upstreamStartPos - max(codingEnd + 1, exon.Start) + 1;

                    if(inCoding)
                    {
                        cc.UpstreamPhase =
                                calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, max(upstreamStartPos, exon.Start));

                        cc.CodingPositionRange[SE_START] = max(max(exon.Start, codingStart), variant.Position);

                        cc.CodingPositionRange[SE_END] = !variant.isInsert() ?
                                min(min(exon.End, codingEnd), variant.EndPosition) : cc.CodingPositionRange[SE_START];

                        // add in any extra coding bases
                        cc.CodingBase = codingBaseCount + cc.CodingPositionRange[SE_START] - max(codingStart, exon.Start) + 1;
                    }
                }
                else if(nextExon != null && variant.altPositionsOverlap(nextExon.Start, nextExon.End))
                {
                    // deal with next exon the next time around
                    continue;
                }
                else
                {
                    cc.RegionType = INTRONIC;

                    int distanceToPrev = upstreamStartPos - exon.End;
                    int distanceToNext = nextExon.Start - upstreamStartPos;

                    if(distanceToPrev < distanceToNext)
                    {
                        cc.ExonRank = exon.Rank;
                        cc.NearestExonDistance = upstreamStartPos - exon.End;
                    }
                    else
                    {
                        cc.ExonRank = nextExon.Rank;
                        cc.NearestExonDistance = upstreamStartPos - nextExon.Start; // -ve for back into previous intron
                    }

                    if(inCoding)
                    {
                        cc.UpstreamPhase = posStrand ? exon.PhaseEnd : nextExon.PhaseEnd;

                        if(distanceToPrev < distanceToNext)
                            cc.CodingBase = codingBaseCount;
                        else
                            cc.CodingBase = codingBaseCount + 1; // first base of next exon
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
                    if(inCoding && positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        codingBaseCount += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

                    prePosExonicBaseCount += exon.baseLength();

                    if(isCodingTrans && exon.Start <= codingStart)
                        postCodingExonicBases += min(codingStart - 1, exon.End) - exon.Start + 1;

                    if(nextExon != null && upstreamStartPos <= nextExon.End)
                        continue;
                }

                boolean withinExon = variant.altPositionsWithin(exon.Start, exon.End);
                if(withinExon || variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    cc.RegionType = EXONIC;
                    cc.ExonRank = exon.Rank;
                    cc.SpansSpliceJunction = !withinExon;

                    prePosExonicBaseCount += exon.End - upstreamStartPos + 1;

                    // bases from end of coding to the position
                    if(isCodingTrans && exon.Start <= codingStart)
                        postCodingExonicBases += min(codingStart - 1, exon.End) - upstreamStartPos + 1;

                    if(inCoding)
                    {
                        cc.UpstreamPhase =
                                calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, min(upstreamStartPos, exon.End));

                        prePosExonicBaseCount += exon.End - upstreamStartPos + 1;

                        if(!variant.isInsert())
                        {
                            cc.CodingPositionRange[SE_START] = max(max(exon.Start, codingStart), variant.Position);
                            cc.CodingPositionRange[SE_END] = min(min(exon.End, codingEnd), variant.EndPosition);

                        }
                        else
                        {
                            cc.CodingPositionRange[SE_END] = min(min(exon.End, codingEnd), variant.EndPosition);
                            cc.CodingPositionRange[SE_START] = cc.CodingPositionRange[SE_END];
                        }

                        // add in any extra coding bases
                        cc.CodingBase = codingBaseCount + min(codingEnd, exon.End) - cc.CodingPositionRange[SE_END] + 1;
                    }
                }
                else if(nextExon != null && variant.altPositionsOverlap(nextExon.Start, nextExon.End))
                {
                    continue;
                }
                else
                {
                    cc.RegionType = INTRONIC;

                    int distanceToPrev = exon.Start - upstreamStartPos;
                    int distanceToNext = upstreamStartPos - nextExon.End;

                    if(distanceToPrev < distanceToNext)
                    {
                        cc.ExonRank = exon.Rank;
                        cc.NearestExonDistance = exon.Start - upstreamStartPos;
                    }
                    else
                    {
                        cc.ExonRank = nextExon.Rank;
                        cc.NearestExonDistance = nextExon.End - upstreamStartPos; // negative
                    }

                    if(inCoding)
                    {
                        cc.UpstreamPhase = posStrand ? exon.PhaseEnd : nextExon.PhaseEnd;

                        if(distanceToPrev < distanceToNext)
                        {
                            cc.CodingBase = codingBaseCount;
                        }
                        else
                        {
                            cc.CodingBase = codingBaseCount + 1; // first base of next exon
                        }
                    }
                }

                break;
            }
        }

        if(cc.CodingType == NON_CODING)
        {
            cc.CodingBase = prePosExonicBaseCount;
        }
        else if(cc.CodingType == UTR_3P)
        {
            cc.CodingBase = postCodingExonicBases;
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

                cc.DeletedCodingBases = ref.length() - alt.length();
                adjustedCodingBases = cc.DeletedCodingBases;
            }

            if(!isCodonMultiple(adjustedCodingBases))
                cc.IsFrameShift = true;
        }
    }

    public static void setPreCodingExonicDistance(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        if(variant.altBasesBelow(transData.CodingStart))
        {
            int codingStart = transData.CodingStart;
            int position = variant.Position; // which ought to be used?

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                if(nextExon != null && variant.altBasesAbove(nextExon.Start - 1))
                    continue;

                if(variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    if(positionWithin(codingStart, exon.Start, exon.End))
                        cc.CodingBase += codingStart - position;
                    else
                        cc.CodingBase += exon.End - position + 1;
                }
                else if(exon.Start > position)
                {
                    if(positionWithin(codingStart, exon.Start, exon.End))
                        cc.CodingBase += codingStart - exon.Start;
                    else
                        cc.CodingBase += exon.baseLength();
                }

                if(positionWithin(codingStart, exon.Start, exon.End))
                    break;
            }
        }
        else
        {
            // similar for position after coding end
            int codingEnd = transData.CodingEnd;
            int position = variant.Position;

            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i >= 1 ? transData.exons().get(i - 1) : null;

                if(nextExon != null && variant.altBasesBelow(nextExon.End + 1))
                    continue;

                if(variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    if(positionWithin(codingEnd, exon.Start, exon.End))
                        cc.CodingBase += position - codingEnd;
                    else
                        cc.CodingBase += position - exon.Start + 1;
                }
                else if(exon.End < position)
                {
                    if(positionWithin(codingEnd, exon.Start, exon.End))
                        cc.CodingBase += exon.End - codingEnd;
                    else
                        cc.CodingBase += exon.baseLength();
                }

                if(positionWithin(codingEnd, exon.Start, exon.End))
                    break;
            }
        }

        cc.CodingBase *= -1;
    }


    private static void setCodingContextOld(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        // initially set the region for this variant and the coding base range and phase if any
        // calculate coding bases prior to the upstream position and add any addition if within an exon
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
                    cc.SpansSpliceJunction = !withinExon;

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
                        cc.NearestExonDistance = upstreamStartPos - exon.End;
                    }
                    else
                    {
                        cc.ExonRank = nextExon.Rank;
                        cc.CodingBase = preExonCodingBases + 1; // first base of next exon
                        cc.NearestExonDistance = upstreamStartPos - nextExon.Start; // -ve for back into previous intron
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
                    cc.SpansSpliceJunction = !withinExon;

                    cc.UpstreamPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, min(upstreamStartPos, exon.End));

                    if(!variant.isInsert())
                    {
                        cc.CodingPositionRange[SE_START] = max(max(exon.Start, codingStart), variant.Position);
                        cc.CodingPositionRange[SE_END] = min(min(exon.End, codingEnd), variant.EndPosition);

                    }
                    else
                    {
                        cc.CodingPositionRange[SE_END] = min(min(exon.End, codingEnd), variant.EndPosition);
                        cc.CodingPositionRange[SE_START] = cc.CodingPositionRange[SE_END];
                    }

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
                        cc.NearestExonDistance = exon.Start - upstreamStartPos;
                    }
                    else
                    {
                        cc.ExonRank = nextExon.Rank;
                        cc.CodingBase = preExonCodingBases + 1; // first base of next exon
                        cc.NearestExonDistance = nextExon.End - upstreamStartPos; // negative
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

                cc.DeletedCodingBases = ref.length() - alt.length();
                adjustedCodingBases = cc.DeletedCodingBases;
            }

            if(!isCodonMultiple(adjustedCodingBases))
                cc.IsFrameShift = true;
        }
    }

    public static CodingContext setPreOrPostCodingContext(final VariantData variant, final TranscriptData transData)
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
                        cc.NearestExonDistance += codingStart - position;
                    else
                        cc.NearestExonDistance += exon.End - position;
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
                        cc.NearestExonDistance += codingStart - exon.Start;
                    else
                        cc.NearestExonDistance += exon.baseLength();
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
                        cc.NearestExonDistance += position - codingEnd;
                    else
                        cc.NearestExonDistance += position - exon.Start;
                }
                else if(nextExon != null && position > nextExon.End && position < exon.Start)
                {
                    cc.ExonRank = transData.posStrand() ? nextExon.Rank : exon.Rank;
                    cc.RegionType = INTRONIC;
                }
                else if(exon.End < position)
                {
                    if(positionWithin(codingEnd, exon.Start, exon.End))
                        cc.NearestExonDistance += exon.End - codingEnd;
                    else
                        cc.NearestExonDistance += exon.baseLength();
                }

                if(positionWithin(codingEnd, exon.Start, exon.End))
                    break;
            }
        }

        // convention is to have negative values if 5'UTR
        if(cc.CodingType == UTR_5P)
        {
            cc.NearestExonDistance *= -1;
        }

        return cc;
    }

    public static CodingContext setNonCodingContext(final int posStart, final int posEnd, final TranscriptData transData)
    {
        CodingContext cc = new CodingContext();

        // can set exon ranks and bases from start of transcript but that's it
        cc.CodingType = NON_CODING;

        // how to define bases for HGVC coding if at all?
        if(transData.posStrand())
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);

                if(posStart > exon.End)
                    continue;

                if(posStart < exon.Start)
                {
                    cc.ExonRank = exon.Rank - 1;
                    cc.RegionType = INTRONIC;

                    cc.NearestExonDistance = posStart - transData.TransStart;

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
            cc.NearestExonDistance = transData.TransEnd - posEnd;

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

    public static void setUpstreamContext(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        cc.RegionType = UPSTREAM;

        // upstream (or downstream but not handled)
        if(transData.posStrand())
        {
            // will be negative
            cc.NearestExonDistance = variant.Position - transData.TransStart;
        }
        else
        {
            cc.NearestExonDistance = transData.TransEnd - variant.Position;
        }
    }

}
