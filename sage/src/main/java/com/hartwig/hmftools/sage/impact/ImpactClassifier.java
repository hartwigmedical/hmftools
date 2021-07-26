package com.hartwig.hmftools.sage.impact;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.START_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.START_CODON;
import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.isStopCodon;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.getCodingBaseRanges;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.FIVE_PRIME_UTR_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.INTRON_VARIANT_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_REGION_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.THREE_PRIME_UTR_EFFECT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRON_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_REGION_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.START_LOST;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_LOST;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SYNONYMOUS_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.UPSTREAM_GENE_VARIANT;
import static com.hartwig.hmftools.sage.impact.CodingContext.determineContext;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.GENE_UPSTREAM_DISTANCE;
import static com.hartwig.hmftools.sage.impact.SpliceClassifier.isWithinSpliceRegion;

import java.util.List;

import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.CodingBaseData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptUtils;
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
        // non-coding transcripts are ignored
        if(transData.CodingStart == null)
            return null;

        int position = variant.Position;

        if(isOutsideTransRange(transData, position))
            return null;

        VariantTransImpact transImpact = checkNonCodingImpact(variant, transData);

        if(transImpact != null)
            return transImpact;

        List<String> consequenceEffects = Lists.newArrayList();

        boolean inSpliceRegion = false;
        int exonRank = 0;

        final List<Integer> nonRefPositions = variant.nonRefPositions();

        // check intron variant
        int exonCount = transData.exons().size();
        for(int i = 0; i < exonCount; ++i)
        {
            ExonData exon = transData.exons().get(i);

            ExonData nextExon = i < exonCount - 1 ? transData.exons().get(i + 1) : null;

            if(!inSpliceRegion && isWithinSpliceRegion(variant, transData, exon))
            {
                inSpliceRegion = true;
                exonRank = exon.Rank;
                addConsequenceEffect(consequenceEffects, mSpliceClassifier.classifyVariant(variant, transData, exon));
            }
            else if(!inSpliceRegion && nextExon != null && isWithinSpliceRegion(variant, transData, nextExon))
            {
                inSpliceRegion = true;
                exonRank = nextExon.Rank;
                addConsequenceEffect(consequenceEffects, mSpliceClassifier.classifyVariant(variant, transData, nextExon));
            }

            if(positionWithin(position, exon.Start, exon.End) || positionWithin(variant.endPosition(), exon.Start, exon.End))
            {
                exonRank = exon.Rank;
                transImpact = classifyExonicPosition(variant, transData, exon);
                break;
            }

            if(nextExon != null && positionsWithin(position, variant.endPosition(), exon.End + 1, nextExon.Start - 1))
            {
                addConsequenceEffect(consequenceEffects, INTRON_VARIANT_EFFECT);
                break;
            }
        }

        if(transImpact == null)
        {
            if(consequenceEffects.isEmpty())
                return null;

            transImpact = new VariantTransImpact(transData, consequenceEffects.get(0));
            consequenceEffects.remove(0);
        }

        for(String effect : consequenceEffects)
        {
            transImpact.addConsequence(effect);
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

    private boolean isOutsideTransRange(final TranscriptData transData, int position)
    {
        if(transData.posStrand())
            return position < transData.TransStart - GENE_UPSTREAM_DISTANCE || position > transData.TransEnd;
        else
            return position < transData.TransStart || position > transData.TransEnd + GENE_UPSTREAM_DISTANCE;
    }

    private VariantTransImpact checkNonCodingImpact(final VariantData variant, final TranscriptData transData)
    {
        int position = variant.Position;

        if(transData.posStrand())
        {
            // check pre-gene region
            if(position >= transData.TransStart - GENE_UPSTREAM_DISTANCE && position < transData.TransStart)
                return new VariantTransImpact(transData, UPSTREAM_GENE_VARIANT);

            // check 5' and 3' UTR
            if(position >= transData.TransStart && position < transData.CodingStart)
                return new VariantTransImpact(transData, FIVE_PRIME_UTR_EFFECT);

            if(position > transData.CodingEnd && position <= transData.TransEnd)
                return new VariantTransImpact(transData, THREE_PRIME_UTR_EFFECT);

        }
        else
        {
            if(position > transData.TransEnd && position <= transData.TransEnd + GENE_UPSTREAM_DISTANCE)
                return new VariantTransImpact(transData, UPSTREAM_GENE_VARIANT);

            if(position >= transData.TransStart && position < transData.CodingStart)
                return new VariantTransImpact(transData, THREE_PRIME_UTR_EFFECT);

            if(position > transData.CodingEnd && position <= transData.TransEnd)
                return new VariantTransImpact(transData, FIVE_PRIME_UTR_EFFECT);
        }

        return null;
    }

    private VariantTransImpact classifyExonicPosition(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        if(variant.isIndel())
        {
            // if only the end of an exon is affected (either strand) then this doesn't change the coding bases
            if(variant.Position == exon.End)
                return null;

            if(variant.isInsert() && variant.endPosition() == exon.Start)
                return null;

            // if the variant crosses the exon boundary then only check the bases which are exonic to decide between frameshift or inframe
            int exonicBases = 0;

            if(variant.isDeletion())
            {
                int firstDelBase = variant.Position + 1;
                int lastDelBase = variant.Position + abs(variant.baseDiff());

                if(abs(exon.Start - variant.Position) < abs(exon.End - variant.Position))
                    exonicBases = lastDelBase - max(firstDelBase, exon.Start) + 1;
                else
                    exonicBases = min(lastDelBase, exon.End) - firstDelBase + 1;
            }
            else
            {
                exonicBases = variant.baseDiff();
            }

            /*
            if(positionWithin(exon.Start, firstDelBase, lastDelBase) && positionWithin(exon.Start + 1, firstDelBase, lastDelBase))
                return new VariantTransImpact(transData, FRAMESHIFT_VARIANT);
            else if(positionWithin(exon.End, firstDelBase, lastDelBase) && positionWithin(exon.End + 1, firstDelBase, lastDelBase))
                return new VariantTransImpact(transData, FRAMESHIFT_VARIANT);
            */

            if((exonicBases % 3) == 0)
            {
                // could use phasing info to set conservative vs disruptive type
                VariantTransImpact transImpact = variant.isDeletion() ?
                        new VariantTransImpact(transData, INFRAME_DELETION) : new VariantTransImpact(transData, INFRAME_INSERTION);

                CodingContext codingContext = determineContext(
                        variant.Chromosome, variant.Position, variant.Ref, variant.Alt, transData, mRefGenome);

                transImpact.setCodingContext(codingContext);

                return transImpact;
            }
            else
            {
                return new VariantTransImpact(transData, FRAMESHIFT_VARIANT);
            }
        }
        else
        {
            CodingContext codingContext = determineContext(
                    variant.Chromosome, variant.Position, variant.Ref, variant.Alt, transData, mRefGenome);

            /*
            // find the start of the codon in which the first base of the variant is in, and likewise the end of the last codon
            int varLength = variant.Ref.length();
            boolean posStrand = transData.posStrand();

            int upstreamStartPos = posStrand ? variant.Position : variant.Position + varLength - 1;
            final CodingBaseData cbData = calcCodingBases(transData, upstreamStartPos);

            // this is 0 if the variant starts on the first base of a codon (phase=1), 1 if it starts on the 2nd, and 2 if on the 3rd/last
            int openCodonBases = getOpenCodonBases(cbData.Phase);

            int codingBaseLen = (int)(Math.ceil((openCodonBases + varLength) / 3.0)) * 3;

            int codonStartPos = posStrand ? upstreamStartPos - openCodonBases : upstreamStartPos + openCodonBases;

            List<int[]> codingRanges = getCodingBaseRanges(transData, codonStartPos, posStrand, codingBaseLen);
            String refCodingBases = mRefGenome.getBaseString(variant.Chromosome, codingRanges);

            String altCodingBases;
            String refAminoAcids;
            String altAminoAcids;

            if(posStrand)
            {
                if(openCodonBases > 0)
                {
                    altCodingBases = refCodingBases.substring(0, openCodonBases) + variant.Alt
                            + refCodingBases.substring(openCodonBases + varLength);
                }
                else
                {
                    altCodingBases = variant.Alt + refCodingBases.substring(varLength);
                }

                refAminoAcids = Codons.aminoAcidFromBases(refCodingBases);
                altAminoAcids = Codons.aminoAcidFromBases(altCodingBases);
            }
            else
            {
                if(openCodonBases > 0)
                {
                    altCodingBases = refCodingBases.substring(0, codingBaseLen - varLength - openCodonBases) + variant.Alt
                            + refCodingBases.substring(codingBaseLen - openCodonBases);
                }
                else
                {
                    altCodingBases = refCodingBases.substring(0, codingBaseLen - varLength) + variant.Alt;
                }

                refAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(refCodingBases));
                altAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(altCodingBases));
            }
            */

            VariantConsequence consequence = codingContext.hasProteinChange() ? MISSENSE_VARIANT : SYNONYMOUS_VARIANT;
            VariantTransImpact transImpact = new VariantTransImpact(transData, consequence);
            transImpact.setCodingContext(codingContext);
            return transImpact;
        }
    }

    private void checkStopStartCodons(final VariantTransImpact transImpact)
    {
        if(!transImpact.hasCodingData())
            return;

        VariantConsequence ssConsequence = checkStopStartCodons(
                transImpact.getCodingContext().WildtypeAA, transImpact.getCodingContext().NovelAA);

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