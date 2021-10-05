package com.hartwig.hmftools.pave;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.START_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.FIVE_PRIME_UTR_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.INTRON_VARIANT_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_REGION_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.THREE_PRIME_UTR_EFFECT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.NON_CODING_TRANSCRIPT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.START_LOST;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_LOST;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SYNONYMOUS_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.UPSTREAM_GENE_VARIANT;
import static com.hartwig.hmftools.pave.PaveConstants.GENE_UPSTREAM_DISTANCE;
import static com.hartwig.hmftools.pave.SpliceClassifier.isWithinSpliceRegion;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.apache.commons.compress.utils.Lists;

public class ImpactClassifier
{
    private final RefGenomeInterface mRefGenome;
    private final SpliceClassifier mSpliceClassifier;

    public ImpactClassifier(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
        mSpliceClassifier = new SpliceClassifier();
    }

    public VariantTransImpact classifyVariant(final VariantData variant, final TranscriptData transData)
    {
        if(isOutsideTransRange(variant, transData))
            return null;

        if(transData.CodingStart == null)
        {
            // check if exonic otherwise ignore
            if(transData.exons().stream().noneMatch(x -> positionsOverlap(variant.Position, variant.EndPosition, x.Start, x.End)))
                return null;
        }

        VariantTransImpact transImpact = new VariantTransImpact(transData);

        CodingContext codingContext = CodingContext.determineContext(variant, transData);
        transImpact.setCodingContext(codingContext);

        if(transData.CodingStart == null)
        {
            transImpact.addConsequence(NON_CODING_TRANSCRIPT_VARIANT);
            return transImpact;
        }

        if(codingContext.isCoding())
        {
            ProteinContext proteinContext = ProteinContext.determineContext(variant, codingContext, transData, mRefGenome);
            transImpact.setProteinContext(proteinContext);
        }

        if(checPrePostCodingImpact(variant, transImpact))
            return transImpact;

        boolean inSpliceRegion = false;
        int exonRank = 0;

        // loop through all exons checking if the variant falls in a splice region and is intronic or exonic
        int exonCount = transData.exons().size();
        for(int i = 0; i < exonCount; ++i)
        {
            ExonData exon = transData.exons().get(i);

            ExonData nextExon = i < exonCount - 1 ? transData.exons().get(i + 1) : null;

            if(!inSpliceRegion && isWithinSpliceRegion(variant, transData, exon))
            {
                inSpliceRegion = true;
                exonRank = exon.Rank;
                transImpact.addConsequence(mSpliceClassifier.classifyVariant(variant, transData, exon));
            }
            else if(!inSpliceRegion && nextExon != null && isWithinSpliceRegion(variant, transData, nextExon))
            {
                inSpliceRegion = true;
                exonRank = nextExon.Rank;
                transImpact.addConsequence(mSpliceClassifier.classifyVariant(variant, transData, nextExon));
            }

            if(positionsWithin(variant.Position, variant.EndPosition, exon.Start, exon.End))
            {
                exonRank = exon.Rank;
                classifyExonicPosition(variant, transImpact, exon);
                break;
            }

            if(nextExon != null && positionsWithin(variant.Position, variant.EndPosition, exon.End + 1, nextExon.Start - 1))
            {
                transImpact.addConsequence(INTRON_VARIANT_EFFECT);
                break;
            }

            // variant which full straddles the splice region
            if(!inSpliceRegion && variant.Position < exon.Start && variant.EndPosition >= exon.Start)
            {
                // TO-DO - confirm that should be splice_donor_variant or splice_acceptor_variant
                transImpact.addConsequence(SPLICE_REGION_EFFECT);
                break;
            }
        }

        if(inSpliceRegion)
        {
            transImpact.markSpliceRegion();
            transImpact.addConsequence(SPLICE_REGION_EFFECT);
        }

        if(exonRank >= 1)
            transImpact.setExonRank(exonRank);

        checkStopStartCodons(transImpact);

        return transImpact;
    }

    private static void addConsequenceEffect(final List<String> consequenceEffects, final String consequenceEffect)
    {
        if(consequenceEffect == null)
            return;

        if(consequenceEffects.contains(consequenceEffect))
            return;

        consequenceEffects.add(consequenceEffect);
    }

    private boolean isOutsideTransRange(final VariantData variant, final TranscriptData transData)
    {
        if(transData.posStrand())
            return variant.EndPosition < transData.TransStart - GENE_UPSTREAM_DISTANCE || variant.Position > transData.TransEnd;
        else
            return variant.EndPosition < transData.TransStart || variant.Position > transData.TransEnd + GENE_UPSTREAM_DISTANCE;
    }

    private boolean checPrePostCodingImpact(final VariantData variant, final VariantTransImpact transImpact)
    {
        final TranscriptData transData = transImpact.TransData;

        if(transData.posStrand())
        {
            // check pre-gene region
            if(variant.Position >= transData.TransStart - GENE_UPSTREAM_DISTANCE && variant.EndPosition < transData.TransStart)
            {
                transImpact.addConsequence(UPSTREAM_GENE_VARIANT);
                return true;
            }

            // check 5' and 3' UTR
            if(variant.Position >= transData.TransStart && variant.EndPosition < transData.CodingStart)
            {
                transImpact.addConsequence(FIVE_PRIME_UTR_EFFECT);
                return true;
            }

            if(variant.Position > transData.CodingEnd && variant.Position <= transData.TransEnd)
            {
                transImpact.addConsequence(THREE_PRIME_UTR_EFFECT);
                return true;
            }
        }
        else
        {
            if(variant.Position > transData.TransEnd && variant.EndPosition <= transData.TransEnd + GENE_UPSTREAM_DISTANCE)
            {
                transImpact.addConsequence(UPSTREAM_GENE_VARIANT);
                return true;
            }

            if(variant.Position >= transData.TransStart && variant.EndPosition < transData.CodingStart)
            {
                transImpact.addConsequence(THREE_PRIME_UTR_EFFECT);
                return true;
            }

            if(variant.Position > transData.CodingEnd && variant.EndPosition <= transData.TransEnd)
            {
                transImpact.addConsequence(FIVE_PRIME_UTR_EFFECT);
                return true;
            }
        }

        return false;
    }

    private void classifyExonicPosition(final VariantData variant, final VariantTransImpact transImpact, final ExonData exon)
    {
        final TranscriptData transData = transImpact.TransData;

        if(variant.isIndel())
        {
            // if only the end of an exon is affected (either strand) then this doesn't change the coding bases
            if(variant.EndPosition == exon.Start || variant.Position == exon.End)
                return;

            // if the variant crosses the exon boundary then only check the bases which are exonic to decide between frameshift or inframe
            int exonicBases = 0;

            if(variant.isDeletion())
            {
                int firstDelBase = variant.nonRefPositions().get(0);
                int lastDelBase = variant.nonRefPositions().get(variant.nonRefPositions().size() - 1);

                if(abs(exon.Start - variant.Position) < abs(exon.End - variant.Position))
                    exonicBases = lastDelBase - max(firstDelBase, exon.Start) + 1;
                else
                    exonicBases = min(lastDelBase, exon.End) - firstDelBase + 1;
            }
            else
            {
                exonicBases = variant.baseDiff();
            }

            if((exonicBases % 3) == 0)
            {
                // could use phasing info to set conservative vs disruptive type

                if(variant.isDeletion())
                    transImpact.addConsequence(INFRAME_DELETION);
                else
                    transImpact.addConsequence(INFRAME_INSERTION);
            }
            else
            {
                transImpact.addConsequence(FRAMESHIFT_VARIANT);
            }
        }
        else
        {
            if(transImpact.getProteinContext().hasProteinChange())
                transImpact.addConsequence(MISSENSE_VARIANT);
            else
                transImpact.addConsequence(SYNONYMOUS_VARIANT);
        }
    }

    private void checkStopStartCodons(final VariantTransImpact transImpact)
    {
        if(!transImpact.hasCodingData())
            return;

        VariantConsequence ssConsequence = checkStopStartCodons(
                transImpact.getProteinContext().WildtypeAA, transImpact.getProteinContext().NovelAA);

        if(ssConsequence != null)
        {
            transImpact.addConsequence(ssConsequence.description());
        }
    }

    public static VariantConsequence checkStopStartCodons(final String refAminoAcids, final String altAminoAcids)
    {
        if(refAminoAcids.isEmpty() || altAminoAcids.isEmpty())
            return null;

        if(refAminoAcids.charAt(0) == START_AMINO_ACID && altAminoAcids.charAt(0) != START_AMINO_ACID)
            return START_LOST;

        if(refAminoAcids.charAt(refAminoAcids.length() - 1) != STOP_AMINO_ACID)
        {
            // has the alt gained a stop codon anywhere
            for(int i = 0; i < altAminoAcids.length(); ++i)
            {
                if(altAminoAcids.charAt(i) == STOP_AMINO_ACID)
                    return STOP_GAINED;
            }
        }
        else
        {
            boolean hasStopCodon = false;
            for(int i = 0; i < altAminoAcids.length(); ++i)
            {
                if(altAminoAcids.charAt(i) == STOP_AMINO_ACID)
                {
                    hasStopCodon = true;
                    break;
                }
            }

            if(!hasStopCodon)
                return STOP_LOST;
        }

        return null;
    }
}