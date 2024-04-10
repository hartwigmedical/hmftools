package com.hartwig.hmftools.pave.impact;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.ENHANCER;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UNKNOWN;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcExonicCodingPhase;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.PaveConstants.PROMOTOR_UPSTREAM_DISTANCE;
import static com.hartwig.hmftools.pave.PaveConstants.PROMOTOR_UPSTREAM_GENE_IDS;
import static com.hartwig.hmftools.pave.impact.ProteinUtils.getOpenCodonBases;
import static com.hartwig.hmftools.pave.impact.SpliceClassifier.isInsertIntoExonStart;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pave.VariantData;

public final class CodingUtils
{
    public static void determineContext(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        cc.Strand = transData.Strand;

        setCodingContext(variant, transData, cc);

        // set pre-coding exonic bases separately since doesn't fit with way exons are iterated through
        if(cc.CodingType == UTR_5P || cc.CodingType == UTR_3P || cc.CodingType == ENHANCER)
            setUtrCodingExonicDistance(variant, transData, cc);

        if(variant.isDeletion() && cc.RegionType == EXONIC)
        {
            // calculate trimmed deleted bases within exons
            String ref = cc.codingRef(variant);
            String alt = cc.codingAlt(variant);
            cc.DeletedCodingBases = ref.length() - alt.length();
        }

        // now coding bases have been clarified, mark any frameshift INDELs
        if(variant.isIndel() && cc.CodingType == CODING && cc.RegionType == EXONIC)
        {
            int adjustedCodingBases;

            if(variant.isInsert())
            {
                adjustedCodingBases = variant.baseDiff();
            }
            else
            {
                adjustedCodingBases = cc.DeletedCodingBases;

                if(cc.NearestExonDistance != 0)
                {
                    // used to set HGVS coding string correctly
                    if(transData.posStrand() && variant.altPositions().contains(transData.CodingEnd))
                        cc.SpansCodingEnd = true;
                    else if(!transData.posStrand() && variant.altPositions().contains(transData.CodingStart))
                        cc.SpansCodingEnd = true;
                }
            }

            if(!isCodonMultiple(adjustedCodingBases))
                cc.IsFrameShift = true;
        }
    }

    private static void setCodingContext(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        // set the following data
        // region type
        // coding type
        // distance to nearest exon (or upstream of transcript start)
        // coding base if within coding or number of exonic bases to coding start or end if 5' or 3' UTR

        boolean inCoding = false;
        boolean isCodingTrans = !transData.nonCoding();
        Integer codingStart = transData.CodingStart;
        Integer codingEnd = transData.CodingEnd;

        if(transData.nonCoding())
            cc.CodingType = NON_CODING;

        if(variant.altBasesBelow(transData.TransStart) || variant.altBasesAbove(transData.TransEnd))
        {
            cc.RegionType = UPSTREAM;

            if(!transData.nonCoding() && PROMOTOR_UPSTREAM_GENE_IDS.contains(transData.GeneId))
                cc.CodingType = ENHANCER;

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

                boolean withinExon = variant.altPositionsWithin(exon.Start, exon.End)
                        || isInsertIntoExonStart(variant, exon, posStrand);

                boolean overlapsExon = variant.altPositionsOverlap(exon.Start, exon.End);

                if(skipExonPos > exon.End)
                {
                    // keep track of coding bases if this exon overlaps the coding region but is before the variant
                    if(inCoding && !withinExon && !overlapsExon && positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        codingBaseCount += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

                    if(variant.Position > exon.End)
                        prePosExonicBaseCount += exon.baseLength();

                    if(nextExon != null && upstreamStartPos >= nextExon.Start)
                        continue;
                }

                // check for fully intronic, fully exonic or spanning

                if(withinExon || overlapsExon)
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

                        // measure distance of a variant into the intron if applicable
                        if(!variant.isInsert())
                        {
                            int exonCodingStart = max(exon.Start, codingStart);
                            int exonCodingEnd = min(exon.End, codingEnd);

                            if(variant.EndPosition > exonCodingEnd)
                                cc.NearestExonDistance = variant.EndPosition - exonCodingEnd;
                            else if(variant.Position < exonCodingStart)
                                cc.NearestExonDistance = variant.Position - exonCodingStart; // -ve back into previous intron

                            cc.SpansCodingStart = variant.Position < codingStart;
                        }
                    }
                    else
                    {
                        // set exonic range and handle variants that span from exon to intron
                        cc.CodingPositionRange[SE_START] = max(exon.Start, variant.Position);

                        cc.CodingPositionRange[SE_END] = !variant.isInsert() ?
                                min(exon.End, variant.EndPosition) : cc.CodingPositionRange[SE_START];

                        if(!variant.isInsert())
                        {
                            if(variant.EndPosition > exon.End)
                                cc.NearestExonDistance = variant.EndPosition - exon.End;
                            else if(variant.Position < exon.Start)
                                cc.NearestExonDistance = variant.Position - exon.Start; // -ve back into previous intron
                        }
                    }
                }
                else if(nextExon != null
                && (variant.altPositionsOverlap(nextExon.Start, nextExon.End)) || isInsertIntoExonStart(variant, nextExon, posStrand))
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
            // negative strand transcripts
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

                boolean withinExon = variant.altPositionsWithin(exon.Start, exon.End)
                        || isInsertIntoExonStart(variant, exon, posStrand);

                boolean overlapsExon = variant.altPositionsOverlap(exon.Start, exon.End);

                if(skipExonPos < exon.Start)
                {
                    if(inCoding && !withinExon && !overlapsExon && positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        codingBaseCount += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

                    if(variant.EndPosition < exon.Start)
                        prePosExonicBaseCount += exon.baseLength();

                    if(nextExon != null && upstreamStartPos <= nextExon.End)
                        continue;
                }

                if(withinExon || overlapsExon)
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

                        if(!variant.isInsert())
                        {
                            int exonCodingStart = max(exon.Start, codingStart);
                            int exonCodingEnd = min(exon.End, codingEnd);

                            if(variant.EndPosition > exonCodingEnd)
                                cc.NearestExonDistance = exonCodingEnd - variant.EndPosition; // -ve back into previous intron
                            else if(variant.Position < exonCodingStart)
                                cc.NearestExonDistance = exonCodingStart - variant.Position;

                            cc.SpansCodingStart = variant.EndPosition > codingEnd;
                        }
                    }
                    else
                    {
                        if(!variant.isInsert())
                        {
                            cc.CodingPositionRange[SE_START] = max(exon.Start, variant.Position);
                            cc.CodingPositionRange[SE_END] = min(exon.End, variant.EndPosition);

                            if(variant.EndPosition > exon.End)
                                cc.NearestExonDistance = exon.End - variant.EndPosition; // -ve back into previous intron
                            else if(variant.Position < exon.Start)
                                cc.NearestExonDistance = exon.Start - variant.Position;
                        }
                        else
                        {
                            cc.CodingPositionRange[SE_END] = min(exon.End, variant.EndPosition);
                            cc.CodingPositionRange[SE_START] = cc.CodingPositionRange[SE_END];
                        }
                    }
                }
                else if(nextExon != null
                && (variant.altPositionsOverlap(nextExon.Start, nextExon.End)) || isInsertIntoExonStart(variant, nextExon, posStrand))
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
    }

    public static void setUtrCodingExonicDistance(final VariantData variant, final TranscriptData transData, final CodingContext cc)
    {
        // for positions in the UTR, determine the exonic distance between the position and the start or end of the coding region
        // 5'UTR position +ve strand transcript, distance from position to coding start
        // 3'UTR position +ve strand, distance from coding end to position
        // 5'UTR position -ve strand, distance from coding end to position
        // 3'UTR position -ve strand, distance from positon to coding start
        boolean posStrand = transData.posStrand();
        boolean posBeforeCodingStart = (posStrand == (cc.CodingType == UTR_5P || cc.CodingType == ENHANCER));

        int codingStart = transData.CodingStart;
        int codingEnd = transData.CodingEnd;
        int position = posStrand ? variant.Position : variant.EndPosition;
        int totalCodingBases = 0;
        boolean codingStartsEndsOnExonBoundary = false;

        for(int i = 0; i < transData.exons().size(); ++i)
        {
            final ExonData exon = transData.exons().get(i);

            if(posBeforeCodingStart)
            {
                if(exon.Start > codingStart)
                    break;

                if(position <= exon.End)
                {
                    // take any portion of exonic bases between the position and coding
                    if(positionWithin(codingStart, exon.Start, exon.End))
                        cc.CodingBase += codingStart - max(position, exon.Start);
                    else
                        cc.CodingBase += exon.End - max(position, exon.Start) + 1;
                }
            }
            else
            {
                if(exon.Start > position)
                    break;

                if(codingEnd <= exon.End)
                {
                    if(positionWithin(codingEnd, exon.Start, exon.End))
                        cc.CodingBase += min(position, exon.End) - codingEnd;
                    else
                        cc.CodingBase += min(position, exon.End) - exon.Start + 1;
                }
            }

            if(positionsOverlap(codingStart, codingEnd, exon.Start, exon.End))
                totalCodingBases += min(codingEnd, exon.End) - max(codingStart, exon.Start) + 1;

            if(codingEnd == exon.End || codingStart == exon.Start)
                codingStartsEndsOnExonBoundary = true;
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
            if(codingStartsEndsOnExonBoundary)
                cc.CodingEndsOnExonBoundary = true;

            if(cc.CodingType == UTR_5P)
                cc.CodingBase = 1; // since zero-based
            else
                cc.CodingBase = totalCodingBases;
        }

        if(cc.CodingType == ENHANCER)
        {
            if(posStrand)
                cc.CodingBase = transData.CodingStart - position;
            else
                cc.CodingBase = position - transData.CodingEnd;

            // exclude this type if too far upstream
            if(cc.CodingBase > PROMOTOR_UPSTREAM_DISTANCE)
            {
                cc.CodingBase = 0;
                cc.CodingType = UNKNOWN;
            }
        }
    }
}
