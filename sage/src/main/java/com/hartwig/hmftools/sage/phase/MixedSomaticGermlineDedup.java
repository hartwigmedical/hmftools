package com.hartwig.hmftools.sage.phase;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcExonicCodingPhase;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;

import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.select.TranscriptSelector;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

public class MixedSomaticGermlineDedup extends BufferedPostProcessor
{
    private static final int MAX_DISTANCE = 10;

    private final TranscriptSelector mSelector;

    public MixedSomaticGermlineDedup(final Consumer<SageVariant> consumer, final List<TranscriptData> transcripts)
    {
        super(MAX_DISTANCE, consumer);
        mSelector = new TranscriptSelector(transcripts);
    }

    @Override
    protected void processSageVariant(final SageVariant variant, final Collection<SageVariant> variants)
    {
        if(variant.isIndel() || !variant.hasLocalPhaseSets())
            return;

        boolean newVariantIsSnv = isPassingSnv(variant);
        boolean newVariantIsMixedGermlineMnv = isMixedGermlineMnv(variant);

        if(!newVariantIsSnv && !newVariantIsMixedGermlineMnv)
            return;

        for(SageVariant other : variants)
        {
            if(other.hasMatchingLps(variant.localPhaseSets()))
            {
                if(newVariantIsSnv && isMixedGermlineMnv(other))
                {
                    process(other, variant);
                }
                else if(newVariantIsMixedGermlineMnv && isPassingSnv(other))
                {
                    process(variant, other);
                }
            }
        }
    }

    private void process(final SageVariant mnv, final SageVariant snv)
    {
        if(longerContainsShorter(snv, mnv))
        {
            snv.mixedGermlineImpact(mnv.mixedGermlineImpact());

            final BaseRegion positionCodon = findCodon(snv.position());

            if(positionCodon != null && keepMnv(positionCodon, snv.variant(), mnv.variant()))
            {
                snv.filters().add(VariantVCF.DEDUP_FILTER);
            }
            else
            {
                mnv.filters().add(VariantVCF.DEDUP_FILTER);
            }
        }
    }

    private BaseRegion findCodon(int position)
    {
        final TranscriptData transcript = mSelector.select(position);

        if(transcript == null || transcript.nonCoding())
            return null;

        if(!positionWithin(position, transcript.CodingStart, transcript.CodingEnd))
            return null;

        // find the codon surrounding this position
        ExonData exon = transcript.exons().stream().filter(x -> positionWithin(position, x.Start, x.End)).findFirst().orElse(null);

        if(exon == null)
            return null;

        int codingPhase = calcExonicCodingPhase(exon, transcript.CodingStart, transcript.CodingEnd, transcript.Strand, position);
        int codonStart = position;
        int codonEnd = position;

        if(codingPhase == PHASE_1)
        {
            codonEnd = min(exon.End, position + 2);
        }
        else if(codingPhase == PHASE_2)
        {
            codonStart = max(exon.Start, position - 1);
            codonEnd = min(exon.End, position + 1);
        }
        else if(codingPhase == PHASE_0)
        {
            codonStart = max(exon.Start, position - 2);
        }

        return new BaseRegion(codonStart, codonEnd);
    }

    private static boolean isPassingSnv(final SageVariant variant)
    {
        return variant.isPassing() && variant.isSnv();
    }

    private static boolean isMixedGermlineMnv(final SageVariant variant)
    {
        return variant.isPassing() && variant.isMnv() && variant.mixedGermlineImpact() > 0;
    }

    static boolean keepMnv(final BaseRegion codon, final VariantHotspot somaticSnv, final VariantHotspot mixedMnv)
    {
        int snvCodonDifferences = codonDifferences(codon, somaticSnv);
        int mnvCodonDifferences = codonDifferences(codon, mixedMnv);

        return mnvCodonDifferences > snvCodonDifferences;
    }

    public static int codonDifferences(final BaseRegion codon, final VariantHotspot variant)
    {
        if(codon.start() > variant.end() || variant.position() > codon.end())
            return 0;

        int overlapStart = max(codon.start(), variant.position());
        int overlapEnd = min(codon.end(), variant.end());

        int difference = 0;
        for(int position = overlapStart; position <= overlapEnd; position++)
        {
            int variantIndex = position - variant.position();

            if(variant.ref().charAt(variantIndex) != variant.alt().charAt(variantIndex))
            {
                difference++;
            }
        }
        return difference;
    }
}
