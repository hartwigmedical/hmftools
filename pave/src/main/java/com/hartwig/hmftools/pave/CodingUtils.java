package com.hartwig.hmftools.pave;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
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
import static com.hartwig.hmftools.pave.ProteinUtils.getOpenCodonBases;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;

public final class CodingUtils
{
    public static void determineContext(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        cc.Strand = transData.Strand;

        setCodingContext(variant, transData, cc);

        // set pre-coding exonic bases separately since doesn't fit with way exons are iterated through
        if(cc.CodingType == UTR_5P || cc.CodingType == UTR_3P)
            setUtrCodingExonicDistance(variant, transData, cc);

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

        if(transData.nonCoding())
            cc.CodingType = NON_CODING;

        if(variant.altBasesBelow(transData.TransStart) || variant.altBasesAbove(transData.TransEnd))
        {
            cc.RegionType = UPSTREAM;
            return;
        }

        boolean posStrand = transData.posStrand();

        if(isCodingTrans)
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
        int codingBaseCount = 0; // accumulated since start of coding region prior to exon closest to variant

        int upstreamStartPos = posStrand ? variant.Position : variant.EndPosition;

        // INDELs don't affect the last base of an exon so need diff criteria to skip it
        int skipExonPos = (!variant.isIndel() == posStrand) ? variant.Position : variant.EndPosition;

        if(posStrand)
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);

                if(i == 0 && isCodingTrans && codingStart == exon.Start)
                {
                    // for rare incomplete transcript where coding start phase is not 1, typically at first base of exon 1
                    int startPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, exon.Start);
                    if(startPhase != PHASE_1)
                        codingBaseCount -= 3 - getOpenCodonBases(startPhase);
                }

                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                if(skipExonPos > exon.End)
                {
                    // keep track of coding bases if this exon overlaps the coding region but is before the variant
                    if(inCoding && positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        codingBaseCount += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

                    prePosExonicBaseCount += exon.baseLength();

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

                if(exon.Rank == 1 && isCodingTrans && codingEnd == exon.End)
                {
                    int startPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, exon.End);
                    if(startPhase != PHASE_1)
                        codingBaseCount -= 3 - getOpenCodonBases(startPhase);
                }

                final ExonData nextExon = i >= 1 ? transData.exons().get(i - 1) : null;

                if(skipExonPos < exon.Start)
                {
                    if(inCoding && positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        codingBaseCount += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

                    prePosExonicBaseCount += exon.baseLength();

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

                    if(inCoding)
                    {
                        cc.UpstreamPhase =
                                calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, min(upstreamStartPos, exon.End));

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

            if(cc.NearestExonDistance < 0)
                cc.CodingBase += 1; // to the next exonic base
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

    public static void setUtrCodingExonicDistance(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        boolean posStrand = transData.posStrand();
        boolean posBeforeCodingStart = (posStrand == (cc.CodingType == UTR_5P));

        int codingStart = transData.CodingStart;
        int codingEnd = transData.CodingEnd;
        int position = posStrand ? variant.Position : variant.EndPosition; // which ought to be used?
        int totalCodingBases = 0;
        boolean codingEndsOnExonBoundary = false;

        for(int i = 0; i < transData.exons().size(); ++i)
        {
            final ExonData exon = transData.exons().get(i);

            if(posBeforeCodingStart)
            {
                // exons before the position
                if(position > exon.End)
                    continue;

                if(exon.Start > codingStart)
                    break;
            }
            else
            {
                if(codingEnd > exon.End)
                    continue;

                if(exon.Start > position)
                    break;
            }

            if(posBeforeCodingStart)
            {
                // take any portion of exonic bases between the position and coding
                if(positionWithin(codingStart, exon.Start, exon.End))
                    cc.CodingBase += codingStart - max(position, exon.Start);
                else
                    cc.CodingBase += exon.End - max(position, exon.Start) + 1;
            }
            else
            {
                if(positionWithin(codingEnd, exon.Start, exon.End))
                    cc.CodingBase += min(position, exon.End) - codingEnd;
                else
                    cc.CodingBase += min(position, exon.End) - exon.Start + 1;
            }

            if(positionsOverlap(codingStart, codingEnd, exon.Start, exon.End))
                totalCodingBases += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

            if(posStrand && codingEnd == exon.End)
                codingEndsOnExonBoundary = true;
            else if(!posStrand && codingStart == exon.Start)
                codingEndsOnExonBoundary = true;
        }

        // push base by 1 if intronic and closest to the next exon
        if(cc.RegionType == INTRONIC)
        {
            // if closer to an exon away from coding start or end, then add a base
            if(cc.CodingType == UTR_5P && cc.NearestExonDistance > 0) // back one base upstream
                ++cc.CodingBase;
            else if(cc.CodingType == UTR_3P && cc.NearestExonDistance < 0) // one base downstream
                ++cc.CodingBase;
        }

        if(cc.CodingBase == 0)
        {
            if(codingEndsOnExonBoundary)
                cc.CodingBase = totalCodingBases;
            else
                cc.CodingBase = 1; // since zero-based
        }
    }
}
