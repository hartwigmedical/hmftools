package com.hartwig.hmftools.sage.impact;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.subtractExact;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_ACCEPTOR_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_DONOR_EFFECT;
import static com.hartwig.hmftools.common.variant.SpliceSites.getAcceptorPosition;
import static com.hartwig.hmftools.common.variant.SpliceSites.getDonorPosition;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_DONOR_VARIANT;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_ACCEPTOR_BASES;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_DONOR_EXON_BASES;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.SPLICE_REGION_DISTANCE;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.VariantConsequence;

public class SpliceClassifier
{
    public SpliceClassifier()
    {
    }

    public static boolean isWithinSpliceRegion(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        int exonStartSpliceRegionStart = exon.Start - SPLICE_REGION_DISTANCE;
        int exonStartSpliceRegionEnd;

        int exonEndSpliceRegionStart;
        int exonEndSpliceRegionEnd = exon.End + SPLICE_REGION_DISTANCE;

        if(transData.Strand == POS_STRAND)
        {
            exonStartSpliceRegionEnd = exon.Start + SPLICE_ACCEPTOR_BASES;
            exonEndSpliceRegionStart = exon.End - SPLICE_DONOR_EXON_BASES;
        }
        else
        {
            exonStartSpliceRegionEnd = exon.Start + SPLICE_DONOR_EXON_BASES;
            exonEndSpliceRegionStart = exon.End - SPLICE_ACCEPTOR_BASES;
        }

        boolean atTransStart = exon.Start == transData.TransStart;
        boolean atTransEnd = exon.End == transData.TransEnd;

        if(!atTransStart && exonStartSpliceRegionStart <= variant.Position && variant.Position <= exonStartSpliceRegionEnd)
            return true;
        else if(!atTransEnd && exonEndSpliceRegionStart <= variant.Position && variant.Position <= exonEndSpliceRegionEnd)
            return true;
        else
            return false;
    }

    public VariantTransImpact classifyVariant(final VariantData variant, final TranscriptData transData, final ExonData exon)
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
            if(variant.isInsert())
            {
                // TODO - for INDELs check repeat sequence vs microhomology condition described above

            }
            else
            {
                // if deletes any coding bases then handle in normal logic
                boolean deletesD5 = false;

                for(int i = 1; i < variant.Ref.length(); ++i) // starts at 1 since the alt exists at this base
                {
                    int position = variant.Position + i;

                    if(position >= exon.Start && position <= exon.End)
                        return null;

                    int donorPos = getDonorPosition(position, donorExonPos, transData.strand());
                    if(donorPos == 5)
                        deletesD5 = true;
                }

                if(!transData.posStrand() && deletesD5)
                {
                    int exonDistanceStart = abs(variant.Position + 1 - donorExonPos);
                    int exonDistanceEnd = abs(variant.Position + variant.Ref.length() - 11 - donorExonPos);
                    int exonDistance = min(exonDistanceStart, exonDistanceEnd);

                    if(exonDistance - variant.Ref.length() <= 5)
                        return new VariantTransImpact(transData, SPLICE_DONOR_EFFECT);
                }
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

                    if(donorPos == -1 || donorPos == 1 || donorPos == 2 || donorPos == 5)
                        return new VariantTransImpact(transData, SPLICE_DONOR_EFFECT);
                }
                else
                {
                    int acceptorPos = getAcceptorPosition(position, acceptorExonPos, transData.strand());

                    if(acceptorPos == 1 || acceptorPos == 2)
                        return new VariantTransImpact(transData, SPLICE_ACCEPTOR_EFFECT);

                    if(acceptorPos == 3)
                    {
                        final char requiredBase = transData.posStrand() ? 'G' : 'C';

                        if(altBase == requiredBase)
                            return new VariantTransImpact(transData, SPLICE_ACCEPTOR_EFFECT);
                    }
                }
            }
        }

        // return new VariantTransImpact(transData, VariantConsequence.SPLICE_REGION_VARIANT);
        return null;
    }
}
