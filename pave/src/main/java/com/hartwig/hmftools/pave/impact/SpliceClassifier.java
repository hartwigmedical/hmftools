package com.hartwig.hmftools.pave.impact;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.variant.SpliceSites.getAcceptorPosition;
import static com.hartwig.hmftools.common.variant.SpliceSites.getDonorPosition;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SPLICE_ACCEPTOR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SPLICE_DONOR;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_ACCEPTOR_END_RANGE;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_ACCEPTOR_POSITIONS;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_DONOR_POSITIONS;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_REGION_EXON_RANGE;
import static com.hartwig.hmftools.pave.PaveConstants.SPLICE_REGION_INTRON_RANGE;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.BASE_CHANGE;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.BASE_SHIFT;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.HOMOLOGY_SHIFT;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.OUTSIDE_RANGE;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.REGION_DELETED;
import static com.hartwig.hmftools.pave.impact.SpliceImpactType.UNKNOWN;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.pave.VariantData;

import org.apache.commons.compress.utils.Lists;

public final class SpliceClassifier
{
    public static boolean isWithinSpliceRegion(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        int exonStartSpliceRegionStart = exon.Start - SPLICE_REGION_INTRON_RANGE;
        int exonStartSpliceRegionEnd = exon.Start + SPLICE_REGION_EXON_RANGE - 1;

        int exonEndSpliceRegionStart = exon.End - SPLICE_REGION_EXON_RANGE + 1;
        int exonEndSpliceRegionEnd = exon.End + SPLICE_REGION_INTRON_RANGE;

        boolean atTransStart = exon.Start == transData.TransStart;
        boolean atTransEnd = exon.End == transData.TransEnd;

        if(!atTransStart)
        {
            if(variant.altPositionsOverlap(exonStartSpliceRegionStart, exonStartSpliceRegionEnd))
                return true;
        }

        if(!atTransEnd)
        {
            if(variant.altPositionsOverlap(exonEndSpliceRegionStart, exonEndSpliceRegionEnd))
                return true;
        }

        return false;
    }

    public static VariantEffect checkStraddlesSpliceRegion(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        if(!variant.isDeletion())
            return null;

        if(exon.Start != transData.TransStart && positionWithin(exon.Start, variant.Position, variant.EndPosition))
        {
            return transData.posStrand() ? SPLICE_ACCEPTOR : SPLICE_DONOR;
        }
        else if(exon.End != transData.TransEnd && positionWithin(exon.End, variant.Position, variant.EndPosition))
        {
            return transData.posStrand() ? SPLICE_DONOR : SPLICE_ACCEPTOR;
        }

        return null;
    }

    public static void classifyVariant(
            final VariantData variant, final VariantTransImpact transImpact, final ExonData exon, final RefGenomeInterface refGenome)
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

        final TranscriptData transData = transImpact.TransData;
        boolean posStrand = transData.posStrand();

        int donorExonPos = posStrand ? exon.End : exon.Start;
        int acceptorExonPos = posStrand ? exon.Start : exon.End;
        boolean isDonorCandidate = abs(variant.Position - donorExonPos) < abs(variant.Position - acceptorExonPos);

        VariantEffect spliceEffect = null;
        SpliceImpactType impactType = UNKNOWN;
        List<Integer> altPositions = Lists.newArrayList();

        if(variant.isIndel())
        {
            // CHECK - for INDELs check repeat sequence vs microhomology?
            impactType = determineIndelSpliceImpact(variant, refGenome, exon, posStrand);

            if(impactType.isDisruptive())
                spliceEffect = isDonorCandidate ? SPLICE_DONOR : SPLICE_ACCEPTOR;

            // gather up affected bases purely for annotation
            if(variant.isDeletion())
            {
                altPositions.addAll(variant.altPositions());
            }
            else
            {
                // won't capture the actual bases affected from the inserted bases
                altPositions.add(variant.Position);
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
                    {
                        spliceEffect = SPLICE_DONOR;
                        altPositions.add(position);
                    }
                }
                else
                {
                    int acceptorPos = getAcceptorPosition(position, acceptorExonPos, transData.strand());

                    if(SPLICE_ACCEPTOR_POSITIONS.contains(acceptorPos))
                    {
                        if(acceptorPos == SPLICE_ACCEPTOR_END_RANGE)
                        {
                            char requiredA3Base = posStrand ? 'G' : 'C';

                            if(altBase == requiredA3Base)
                            {
                                spliceEffect = SPLICE_ACCEPTOR;
                                altPositions.add(position);
                            }
                        }
                        else
                        {
                            spliceEffect = SPLICE_ACCEPTOR;
                            altPositions.add(position);
                        }
                    }
                }
            }

            if(spliceEffect != null)
                impactType = BASE_CHANGE;
        }

        if(spliceEffect != null)
        {
            transImpact.addEffect(spliceEffect);
            transImpact.setSpliceImpactType(impactType);

            // record the splice bases affected by this variant
            StringJoiner baseInfo = new StringJoiner(ITEM_DELIM);

            for(Integer position : altPositions)
            {
                if(isDonorCandidate)
                {
                    int donorPos = getDonorPosition(position, donorExonPos, transData.strand());

                    if(SPLICE_DONOR_POSITIONS.contains(donorPos))
                        baseInfo.add(String.format("D%d", donorPos));
                }
                else
                {
                    int acceptorPos = getAcceptorPosition(position, acceptorExonPos, transData.strand());

                    if(SPLICE_ACCEPTOR_POSITIONS.contains(acceptorPos))
                        baseInfo.add(String.format("A%d", acceptorPos));
                }
            }

            transImpact.codingContext().SpliceDonorAcceptorBases = baseInfo.toString();
        }
    }

    public static SpliceImpactType determineIndelSpliceImpact(
            final VariantData variant, final RefGenomeInterface refGenome, final ExonData exon, boolean posStrand)
    {
        if(variant.isBaseChange())
            return BASE_CHANGE;

        int donorExonPos = posStrand ? exon.End : exon.Start;
        int acceptorExonPos = posStrand ? exon.Start : exon.End;
        boolean isDonorCandidate = abs(variant.Position - donorExonPos) < abs(variant.Position - acceptorExonPos);
        int exonBoundary = (isDonorCandidate == posStrand) ? exon.End : exon.Start;

        // work out offset for A1-3 and D-1 to D5 relative to the exon boundary
        // in this check 0 = last base of exon, -1 the first of the intron and 5 is the 5th intronic base
        int rangeStart;
        int rangeEnd;

        if(posStrand)
        {
            rangeStart = isDonorCandidate ? 0 : -3;
            rangeEnd = isDonorCandidate ? 5 : -1;
        }
        else
        {
            rangeStart = isDonorCandidate ? -5 : 1;
            rangeEnd = isDonorCandidate ? 0 : 3;
        }

        int posRangeStart = exonBoundary + rangeStart;
        int posRangeEnd = exonBoundary + rangeEnd;

        String refSpliceBases = refGenome.getBaseString(variant.Chromosome, posRangeStart, posRangeEnd);

        // form the alt bases for this same splice region
        String altSpliceBases = "";

        if(variant.isInsert())
        {
            // should have been checked already that the variant does impact the splice regions but check again
            if(posStrand)
            {
                // +ve strand: exon end / splice donor is 30, D-1 to D5 is 30-35, insert must be 30-34 to impact
                // +ve strand: exon start / splice acceptor is 20, A3-A1 is 17-19, insert must be 17-18 to impact
                if(variant.Position >= posRangeEnd || variant.EndPosition <= posRangeStart)
                    return OUTSIDE_RANGE;
            }
            else
            {
                // -ve strand: exon end / splice donor is 20, D-1 to D5 is 15-20, insert must be 15-19 to impact
                // -ve strand: exon start / splice acceptor is 30, A3-A1 is 31-33, insert must be 30 or more to impact
                if(variant.Position >= posRangeEnd || variant.EndPosition <= posRangeStart)
                    return OUTSIDE_RANGE;
            }

            // ignore inserting bases into the acceptor create matching homology - let the realigned variant dictate the net impact
            if(!isDonorCandidate)
                return BASE_SHIFT;

            int preInsertBases = variant.Position - posRangeStart;

            if(preInsertBases > 0)
                altSpliceBases = refGenome.getBaseString(variant.Chromosome, posRangeStart, posRangeStart + preInsertBases - 1);

            // special case of last exon base before acceptor on neg strand
            if(!posStrand && variant.Position == exonBoundary && !isDonorCandidate)
                altSpliceBases += variant.Alt.substring(1);
            else
                altSpliceBases += variant.Alt;

            // eg say insert ref is at 0 (ie last exon base) and inserts 2 bases then need the next 3 ref bases
            int nextRefBase = variant.Position + variant.baseDiff() + 1;
            int postInsertBases = posRangeEnd - nextRefBase;

            if(postInsertBases >= 0)
                altSpliceBases += refGenome.getBaseString(variant.Chromosome, variant.EndPosition, variant.EndPosition + postInsertBases);

            // if too many bases have been retrieved, take the last X to match the length of the ref
            if(altSpliceBases.length() > refSpliceBases.length())
            {
                if(isDonorCandidate == posStrand)
                    altSpliceBases = altSpliceBases.substring(0, refSpliceBases.length());
                else
                    altSpliceBases = altSpliceBases.substring(altSpliceBases.length() - refSpliceBases.length());
            }
        }
        else
        {
            // check for no impact
            if(variant.Position >= posRangeEnd || variant.EndPosition <= posRangeStart)
                return OUTSIDE_RANGE;

            // check for the entire region being deleted
            if(variant.Position < posRangeStart && variant.EndPosition > posRangeEnd)
                return REGION_DELETED;

            // ignore pulling upstream exonic bases to create matching homology - let the realigned variant dictate the net impact
            if(variant.Position < exonBoundary)
                return BASE_SHIFT;

            // otherwise the DEL partially overlaps the splice region and so ref bases need only be pulled from one side or the other
            if(variant.Position >= posRangeStart)
            {
                // start of region is preserved, so take the ref portion up to the first deleted base
                int preDelBases = variant.Position - posRangeStart;

                if(preDelBases >= 0)
                    altSpliceBases = refGenome.getBaseString(variant.Chromosome, variant.Position - preDelBases, variant.Position);

                int postDelBases = posRangeEnd - variant.Position;

                if(postDelBases > 0)
                    altSpliceBases += refGenome.getBaseString(variant.Chromosome, variant.EndPosition, variant.EndPosition + postDelBases - 1);
            }
            else
            {
                // a realigned variant cannot then consider how exonic bases can be used to creating matching homology
                if(variant.isRealignedVariant())
                    return BASE_CHANGE;

                // the DEL starts before the splice range and ends in it, so the comparison alt bases are all pulled from lower positions only
                int preDelBases = variant.EndPosition - posRangeStart;
                altSpliceBases = refGenome.getBaseString(variant.Chromosome, variant.Position - (preDelBases - 1), variant.Position);

                int postDelBases = posRangeEnd - variant.EndPosition;

                if(postDelBases >= 0)
                    altSpliceBases += refGenome.getBaseString(variant.Chromosome, variant.EndPosition, variant.EndPosition + postDelBases);
            }
        }

        if(altSpliceBases.length() != refSpliceBases.length())
        {
            PV_LOGGER.error("splice base mismatch: ref({}) vs alt({})", refSpliceBases, altSpliceBases);
            return UNKNOWN;
        }

        if(!isDonorCandidate)
        {
            // the A3 base of the acceptor region is only check if it becomes a G (or C on -ve strand)
            char refA3Base = posStrand ? refSpliceBases.charAt(0) : refSpliceBases.charAt(2);
            char altA3Base = posStrand ? altSpliceBases.charAt(0) : altSpliceBases.charAt(2);
            char requiredA3Base = posStrand ? 'G' : 'C';

            if(refA3Base == altA3Base || altA3Base != requiredA3Base)
            {
                // exclude this base from comparison
                refSpliceBases = posStrand ? refSpliceBases.substring(1, 3) : refSpliceBases.substring(0, 2);
                altSpliceBases = posStrand ? altSpliceBases.substring(1, 3) : altSpliceBases.substring(0, 2);

                if(refSpliceBases.equals(altSpliceBases))
                    return OUTSIDE_RANGE;
            }
        }

        return refSpliceBases.equals(altSpliceBases) ? HOMOLOGY_SHIFT : BASE_SHIFT;
    }

    public static boolean isInsertIntoExonStart(final VariantData variant, final ExonData exon, boolean posStrand)
    {
        if(!variant.isInsert())
            return false;

        // check if the inserted bases will be incorporated into the coding sequence
        if(posStrand)
            return variant.Position == exon.Start - 1;
        else
            return variant.Position == exon.End;
    }

}
