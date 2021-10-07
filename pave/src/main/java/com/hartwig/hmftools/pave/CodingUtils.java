package com.hartwig.hmftools.pave;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcExonicCodingPhase;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.NonCodingContext.determineNonCodingContext;
import static com.hartwig.hmftools.pave.NonCodingContext.determinePreOrPostCodingContext;
import static com.hartwig.hmftools.pave.NonCodingContext.determineUpstreamContext;

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

        if(transData.posStrand())
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                if(upstreamStartPos > exon.End)
                {
                    // keep track of coding bases if this exon overlaps the coding region but is before the variant
                    if(positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        preExonCodingBases += min(codingEnd, exon.End) - max(codingStart, exon.Start);

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

                    cc.UpstreamPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, upstreamStartPos);

                    // add in any extra coding bases
                    cc.CodingBase = preExonCodingBases + upstreamStartPos - max(codingStart, exon.Start) + 1;

                    cc.CodingPositionRange[SE_START] = max(max(exon.Start, codingStart), variant.Position);
                    cc.CodingPositionRange[SE_END] = min(min(exon.End, codingEnd), variant.EndPosition);

                    // record deleted coding bases since harder to reconstruct later on
                    if(variant.isDeletion())
                    {
                        int firstDelBase = variant.Position + 1;
                        int lastDelBase = variant.EndPosition - 1;
                        int delCodingBaseStart = max(max(exon.Start, codingStart), firstDelBase);
                        int delCodingBaseEnd = min(min(exon.End, codingEnd), lastDelBase);

                        cc.DeletedCodingBases = delCodingBaseEnd - delCodingBaseStart + 1;
                    }
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

                if(upstreamStartPos < exon.Start)
                {
                    if(positionsOverlap(exon.Start, exon.End, codingStart, codingEnd))
                        preExonCodingBases += min(codingEnd, exon.End) - max(codingStart, exon.Start);

                    if(nextExon != null && upstreamStartPos <= nextExon.End)
                        continue;
                }

                boolean withinExon = variant.altPositionsWithin(exon.Start, exon.End);
                if(withinExon || variant.altPositionsOverlap(exon.Start, exon.End))
                {
                    cc.RegionType = EXONIC;
                    cc.ExonRank = exon.Rank;
                    cc.SpansSpiceJunction = !withinExon;

                    cc.UpstreamPhase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, upstreamStartPos);

                    // add in any extra coding bases
                    cc.CodingBase = preExonCodingBases + min(codingEnd, exon.End) - upstreamStartPos + 1;

                    cc.CodingPositionRange[SE_START] = max(max(exon.Start, codingStart), variant.Position);
                    cc.CodingPositionRange[SE_END] = min(min(exon.End, codingEnd), variant.EndPosition);

                    // record deleted coding bases since harder to reconstruct later on
                    if(variant.isDeletion())
                    {
                        int firstDelBase = variant.Position + 1;
                        int lastDelBase = variant.EndPosition - 1;
                        int delCodingBaseStart = max(max(exon.Start, codingStart), firstDelBase);
                        int delCodingBaseEnd = min(min(exon.End, codingEnd), lastDelBase);

                        cc.DeletedCodingBases = delCodingBaseEnd - delCodingBaseStart + 1;
                    }
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

        return cc;
    }

    public static String getExtraBases(
            final TranscriptData transData, final RefGenomeInterface refGenome, final String chromosome,
            final ExonData currentExon, int startPos, int requiredBases, boolean searchUp)
    {
        String extraBases = "";

        if(searchUp)
        {
            int currentExonBases = min(requiredBases, currentExon.End - startPos);

            if(currentExonBases > 0)
            {
                extraBases = refGenome.getBaseString(chromosome, startPos + 1, startPos + 1 + (currentExonBases - 1));
            }

            requiredBases -= currentExonBases;

            if(requiredBases > 0)
            {
                int nextExonRank = transData.posStrand() ? currentExon.Rank + 1 : currentExon.Rank - 1;

                ExonData nextExon = transData.exons().stream().filter(x -> x.Rank == nextExonRank).findFirst().orElse(null);

                String nextExonBases = refGenome.getBaseString(chromosome, nextExon.Start, nextExon.Start + (requiredBases - 1));
                extraBases += nextExonBases;
            }
        }
        else
        {
            // search in lower positions
            int currentExonBases = min(requiredBases, startPos - currentExon.Start);

            if(currentExonBases > 0)
            {
                extraBases = refGenome.getBaseString(chromosome, startPos - currentExonBases, startPos - 1);
            }

            requiredBases -= currentExonBases;

            if(requiredBases > 0)
            {
                int nextExonRank = transData.posStrand() ? currentExon.Rank - 1 : currentExon.Rank + 1;

                ExonData nextExon = transData.exons().stream().filter(x -> x.Rank == nextExonRank).findFirst().orElse(null);

                String nextExonBases = refGenome.getBaseString(chromosome, nextExon.End - (requiredBases - 1), nextExon.End);
                extraBases = nextExonBases + extraBases;
            }
        }

        return extraBases;
    }

}
