package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.codingBaseLength;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;

public final class NonCodingContext
{
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
