package com.hartwig.hmftools.pave;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.START_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FIVE_PRIME_UTR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.INTRONIC;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.NON_CODING_TRANSCRIPT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.THREE_PRIME_UTR;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.UPSTREAM_GENE;
import static com.hartwig.hmftools.pave.PaveConstants.GENE_UPSTREAM_DISTANCE;
import static com.hartwig.hmftools.pave.SpliceClassifier.checkStraddlesSpliceRegion;
import static com.hartwig.hmftools.pave.SpliceClassifier.isWithinSpliceRegion;

import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

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
            transImpact.addEffect(NON_CODING_TRANSCRIPT);
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
                transImpact.addEffect(mSpliceClassifier.classifyVariant(variant, transData, exon));
            }
            else if(!inSpliceRegion && nextExon != null && isWithinSpliceRegion(variant, transData, nextExon))
            {
                inSpliceRegion = true;
                exonRank = nextExon.Rank;
                transImpact.addEffect(mSpliceClassifier.classifyVariant(variant, transData, nextExon));
            }

            if(positionsWithin(variant.Position, variant.EndPosition, exon.Start, exon.End))
            {
                exonRank = exon.Rank;
                classifyExonicPosition(variant, transImpact, exon);
                break;
            }

            if(nextExon != null && positionsWithin(variant.Position, variant.EndPosition, exon.End + 1, nextExon.Start - 1))
            {
                transImpact.addEffect(INTRONIC);
                break;
            }

            // variant which full straddles the splice region
            if(!inSpliceRegion)
            {
                VariantEffect effect = checkStraddlesSpliceRegion(variant, transData, exon);
                if(effect != null)
                {
                    transImpact.addEffect(effect);
                    break;
                }
            }
        }

        if(inSpliceRegion)
        {
            transImpact.markSpliceRegion();
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
                transImpact.addEffect(UPSTREAM_GENE);
                return true;
            }

            // check 5' and 3' UTR
            if(variant.Position >= transData.TransStart && variant.EndPosition < transData.CodingStart)
            {
                transImpact.addEffect(FIVE_PRIME_UTR);
                return true;
            }

            if(variant.Position > transData.CodingEnd && variant.Position <= transData.TransEnd)
            {
                transImpact.addEffect(THREE_PRIME_UTR);
                return true;
            }
        }
        else
        {
            if(variant.Position > transData.TransEnd && variant.EndPosition <= transData.TransEnd + GENE_UPSTREAM_DISTANCE)
            {
                transImpact.addEffect(UPSTREAM_GENE);
                return true;
            }

            if(variant.Position >= transData.TransStart && variant.EndPosition < transData.CodingStart)
            {
                transImpact.addEffect(THREE_PRIME_UTR);
                return true;
            }

            if(variant.Position > transData.CodingEnd && variant.EndPosition <= transData.TransEnd)
            {
                transImpact.addEffect(FIVE_PRIME_UTR);
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
                    transImpact.addEffect(INFRAME_DELETION);
                else
                    transImpact.addEffect(INFRAME_INSERTION);
            }
            else
            {
                transImpact.addEffect(FRAMESHIFT);
            }
        }
        else
        {
            if(transImpact.hasCodingData() && transImpact.getProteinContext().hasProteinChange())
                transImpact.addEffect(MISSENSE);
            else
                transImpact.addEffect(SYNONYMOUS);
        }
    }

    private void checkStopStartCodons(final VariantTransImpact transImpact)
    {
        if(!transImpact.hasCodingData())
            return;

        VariantEffect ssEffect = checkStopStartCodons(
                transImpact.getProteinContext().WildtypeAA, transImpact.getProteinContext().NovelAA);

        if(ssEffect != null)
        {
            transImpact.addEffect(ssEffect);
        }
    }

    public static VariantEffect checkStopStartCodons(final String refAminoAcids, final String altAminoAcids)
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