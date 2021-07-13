package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRON_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.UPSTREAM_GENE_VARIANT;
import static com.hartwig.hmftools.sage.impact.ImpactConstants.GENE_UPSTREAM_DISTANCE;
import static com.hartwig.hmftools.sage.impact.SpliceClassifier.isWithinSpliceRegion;

import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class ImpactClassifier
{
    private final RefGenomeInterface mRefGenome;
    private final SpliceClassifier mSpliceClassifier;

    public ImpactClassifier(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
        mSpliceClassifier = new SpliceClassifier(mRefGenome);
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

        // check intron variant
        for(int i = 0; i < transData.exons().size() - 1; ++i)
        {
            ExonData exon = transData.exons().get(i);
            ExonData nextExon = transData.exons().get(i + 1);

            if(isWithinSpliceRegion(variant, transData, exon) || isWithinSpliceRegion(variant, transData, nextExon))
            {
                return mSpliceClassifier.classifyVariant(variant, transData, exon);
            }

            if(exon.Start <= position && position <= exon.End)
            {
                return classifyExonicPosition(variant, transData, exon);
            }

            if(position > exon.End && position < nextExon.Start)
            {
                return new VariantTransImpact(transData, INTRON_VARIANT);
            }
        }

        return null;
    }

    private boolean isOutsideTransRange(final TranscriptData transData, int position)
    {
        if(transData.Strand == POS_STRAND)
            return position < transData.TransStart - GENE_UPSTREAM_DISTANCE || position > transData.TransEnd;
        else
            return position < transData.TransStart || position > transData.TransEnd + GENE_UPSTREAM_DISTANCE;
    }

    private VariantTransImpact checkNonCodingImpact(final VariantData variant, final TranscriptData transData)
    {
        int position = variant.Position;

        if(transData.Strand == POS_STRAND)
        {
            // check pre-gene region
            if(position >= transData.TransStart - GENE_UPSTREAM_DISTANCE && position < transData.TransStart)
                return new VariantTransImpact(transData, UPSTREAM_GENE_VARIANT);

            // check 5' and 3' UTR
            if(position >= transData.TransStart && position < transData.CodingStart)
                return new VariantTransImpact(transData, "5_prime_UTR_variant");

            if(position > transData.CodingEnd && position <= transData.TransEnd)
                return new VariantTransImpact(transData, "3_prime_UTR_variant");

        }
        else
        {
            if(position > transData.TransEnd && position <= transData.TransEnd)
                return new VariantTransImpact(transData, UPSTREAM_GENE_VARIANT);

            if(position >= transData.TransStart && position < transData.CodingStart)
                return new VariantTransImpact(transData, "3_prime_UTR_variant");

            if(position > transData.CodingEnd && position <= transData.TransEnd)
                return new VariantTransImpact(transData, "5_prime_UTR_variant");
        }

        return null;
    }

    private VariantTransImpact classifyExonicPosition(final VariantData variant, final TranscriptData transData, final ExonData exon)
    {
        // TEMP:
        return new VariantTransImpact(transData, MISSENSE_VARIANT);

        // return null;
    }
}