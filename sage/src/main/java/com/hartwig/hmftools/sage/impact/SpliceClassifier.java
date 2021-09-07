package com.hartwig.hmftools.sage.impact;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_ACCEPTOR_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_DONOR_EFFECT;
import static com.hartwig.hmftools.common.variant.SpliceSites.getAcceptorPosition;
import static com.hartwig.hmftools.common.variant.SpliceSites.getDonorPosition;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_ACCEPTOR_END_RANGE;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_ACCEPTOR_START_RANGE;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_DONOR_END_RANGE;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_DONOR_START_RANGE;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_REGION_EXON_RANGE;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_REGION_INTRON_RANGE;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;

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

        List<Integer> nonRefPositions = variant.nonRefPositions();

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

        return false;
    }

    public String classifyVariant(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        /* Rules:
            - mark any variant hitting slice donor sites (D-1,D+1,D+2,D+5) as splice_donor_variant
            - mark any variant hitting slice acceptor sites (A+1;A+2; A+3 if ALT=G only, ALT=C neg-strand) as splice_acceptor_variant
            - INDELS at D+2 (+ve strand only), D+3 and D+4, D+5 (-ve strand only) should be treated as splice_donor_variant
                UNLESS the repeat sequence = microhomology and the 5th base is not altered
            - DELs near DONOR sites on the negative strand may also delete the D+5 position if exonDistance-length(REF)<=5 and
                should also be treated as splice_donor_variant
            - SNV may also be synonymous, missense or nonsense if at D-1 (ranked by NONSENSE,SPLICE,MISSENSE,SYNONYMOUS
        */

        int donorExonPos = transData.posStrand() ? exon.End : exon.Start;
        int acceptorExonPos = transData.posStrand() ? exon.Start : exon.End;
        boolean isDonorCandidate = abs(variant.Position - donorExonPos) < abs(variant.Position - acceptorExonPos);

        if(variant.isIndel())
        {
            List<Integer> nonRefPositions = variant.nonRefPositions();

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

                    if(donorPos >= SPLICE_DONOR_START_RANGE && donorPos <= SPLICE_DONOR_END_RANGE)
                        return SPLICE_DONOR_EFFECT;

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

                    //if(donorPos == -1 || donorPos == 1 || donorPos == 2 || donorPos == 5)

                    if(donorPos >= SPLICE_DONOR_START_RANGE && donorPos <= SPLICE_DONOR_END_RANGE)
                        return SPLICE_DONOR_EFFECT;
                }
                else
                {
                    int acceptorPos = getAcceptorPosition(position, acceptorExonPos, transData.strand());

                    if(acceptorPos >= SPLICE_ACCEPTOR_START_RANGE && acceptorPos <= SPLICE_ACCEPTOR_END_RANGE)
                    {
                        if(acceptorPos == SPLICE_ACCEPTOR_END_RANGE)
                        {
                            final char requiredBase = transData.posStrand() ? 'G' : 'C';

                            if(altBase == requiredBase)
                                return SPLICE_ACCEPTOR_EFFECT;
                        }
                        else
                        {
                            return SPLICE_ACCEPTOR_EFFECT;
                        }
                    }
                }
            }
        }

        return null;
    }
}
