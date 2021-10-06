package com.hartwig.hmftools.pave;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.variant.SpliceSites.getAcceptorPosition;
import static com.hartwig.hmftools.common.variant.SpliceSites.getDonorPosition;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SPLICE_ACCEPTOR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SPLICE_DONOR;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_ACCEPTOR_END_RANGE;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_ACCEPTOR_POSITIONS;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_DONOR_POSITIONS;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_REGION_EXON_RANGE;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_REGION_INTRON_RANGE;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

public class SpliceClassifier
{
    public SpliceClassifier() {}

    public static boolean isWithinSpliceRegion(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        int exonStartSpliceRegionStart = exon.Start - SPLICE_REGION_INTRON_RANGE;
        int exonStartSpliceRegionEnd = exon.Start + SPLICE_REGION_EXON_RANGE - 1;

        int exonEndSpliceRegionStart = exon.End - SPLICE_REGION_EXON_RANGE + 1;
        int exonEndSpliceRegionEnd = exon.End + SPLICE_REGION_INTRON_RANGE;

        boolean atTransStart = exon.Start == transData.TransStart;
        boolean atTransEnd = exon.End == transData.TransEnd;

        List<Integer> nonRefPositions = variant.altPositions();

        /* doesn't handle INS correctly
        if(!atTransStart)
        {
            if(nonRefPositions.stream().anyMatch(x -> positionWithin(x, exonStartSpliceRegionStart, exonStartSpliceRegionEnd)))
                return true;
        }

        if(!atTransEnd)
        {
            if(nonRefPositions.stream().anyMatch(x -> positionWithin(x, exonEndSpliceRegionStart, exonEndSpliceRegionEnd)))
                return true;
        }
        */

        if(!atTransStart)
        {
            if(positionsOverlap(variant.Position, variant.EndPosition, exonStartSpliceRegionStart, exonStartSpliceRegionEnd))
                return true;
        }

        if(!atTransEnd)
        {
            if(positionsOverlap(variant.Position, variant.EndPosition, exonEndSpliceRegionStart, exonEndSpliceRegionEnd))
                return true;
        }

        return false;
    }

    public static VariantEffect checkStraddlesSpliceRegion(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        if(positionWithin(exon.Start, variant.Position, variant.EndPosition))
        {
            return transData.posStrand() ? SPLICE_ACCEPTOR : SPLICE_DONOR;
        }
        else if(positionWithin(exon.End, variant.Position, variant.EndPosition))
        {
            return transData.posStrand() ? SPLICE_DONOR : SPLICE_ACCEPTOR;
        }

        return null;
    }

    public VariantEffect classifyVariant(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        /* Rules:
            - mark any variant hitting splice donor sites (D-1,D+1,D+2,D+5) as splice_donor_variant
            - mark any variant hitting splice acceptor sites (A+1;A+2; A+3 if ALT=G only, ALT=C neg-strand) as splice_acceptor_variant
            - INDELS at D+2 (+ve strand only), D+3 and D+4, D+5 (-ve strand only) should be treated as splice_donor_variant
                UNLESS the repeat sequence = microhomology and the 5th base is not altered
            - DELs near DONOR sites on the negative strand may also delete the D+5 position if exonDistance-length(REF)<=5 and
                should also be treated as splice_donor_variant
            - SNV may also be synonymous, missense or nonsense if at D-1 (ranked by NONSENSE,SPLICE,MISSENSE,SYNONYMOUS)
        */

        int donorExonPos = transData.posStrand() ? exon.End : exon.Start;
        int acceptorExonPos = transData.posStrand() ? exon.Start : exon.End;
        boolean isDonorCandidate = abs(variant.Position - donorExonPos) < abs(variant.Position - acceptorExonPos);

        if(variant.isIndel())
        {
            List<Integer> nonRefPositions = variant.altPositions();

            if(variant.isInsert())
            {
                // TODO - for INDELs check repeat sequence vs microhomology condition described above

            }
            else
            {
                // if deletes any coding bases then handle in normal logic
                //boolean deletesD5 = false;

                for(Integer position : nonRefPositions)
                {
                    int donorPos = getDonorPosition(position, donorExonPos, transData.strand());

                    if(SPLICE_DONOR_POSITIONS.contains(donorPos))
                        return SPLICE_DONOR;

                    //if(donorPos == SPLICE_DONOR_END_RANGE)
                    //    deletesD5 = true;
                }

                /*
                if(!transData.posStrand() && deletesD5)
                {
                    int exonDistanceStart = abs(variant.Position + 1 - donorExonPos);
                    int exonDistanceEnd = abs(variant.Position + variant.Ref.length() - 11 - donorExonPos);
                    int exonDistance = min(exonDistanceStart, exonDistanceEnd);

                    if(exonDistance - variant.Ref.length() <= SPLICE_DONOR_END_RANGE)
                        return new VariantTransImpact(transData, SPLICE_DONOR_EFFECT);
                }
                */
            }
        }
        else
        {
            for(int i = 0; i < variant.Ref.length(); ++i)
            {
                int position = variant.Position + i;
                char altBase = variant.Alt.charAt(i);

                if(isDonorCandidate)
                {
                    int donorPos = getDonorPosition(position, donorExonPos, transData.strand());

                    if(SPLICE_DONOR_POSITIONS.contains(donorPos))
                        return SPLICE_DONOR;
                }
                else
                {
                    int acceptorPos = getAcceptorPosition(position, acceptorExonPos, transData.strand());

                    if(SPLICE_ACCEPTOR_POSITIONS.contains(acceptorPos))
                    {
                        if(acceptorPos == SPLICE_ACCEPTOR_END_RANGE)
                        {
                            final char requiredBase = transData.posStrand() ? 'G' : 'C';

                            if(altBase == requiredBase)
                                return SPLICE_ACCEPTOR;
                        }
                        else
                        {
                            return SPLICE_ACCEPTOR;
                        }
                    }
                }
            }
        }

        return null;
    }
}
