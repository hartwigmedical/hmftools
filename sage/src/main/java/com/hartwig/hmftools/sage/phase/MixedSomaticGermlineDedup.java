package com.hartwig.hmftools.sage.phase;

import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.select.TranscriptRegionSelector;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import org.jetbrains.annotations.NotNull;

public class MixedSomaticGermlineDedup extends BufferedPostProcessor
{
    private static final int MAX_DISTANCE = 10;
    
    private final TranscriptRegionSelector mSelector;

    public MixedSomaticGermlineDedup(final Consumer<SageVariant> consumer, final List<HmfTranscriptRegion> transcripts)
    {
        super(MAX_DISTANCE, consumer);
        mSelector = new TranscriptRegionSelector(transcripts);
    }

    @Override
    protected void processSageVariant(@NotNull final SageVariant newVariant, @NotNull final Collection<SageVariant> buffer)
    {
        int lps = newVariant.localPhaseSet();

        if(!newVariant.isIndel() && lps > 0)
        {

            boolean newVariantIsSnv = isPassingSnv(newVariant);
            boolean newVariantIsMixedGermlineMnv = isMixedGermlineMnv(newVariant);

            if(newVariantIsSnv || newVariantIsMixedGermlineMnv)
            {
                for(SageVariant oldVariant : buffer)
                {
                    if(oldVariant.localPhaseSet() == lps)
                    {
                        if(newVariantIsSnv && isMixedGermlineMnv(oldVariant))
                        {
                            process(oldVariant, newVariant);
                        }
                        else if(newVariantIsMixedGermlineMnv && isPassingSnv(oldVariant))
                        {
                            process(newVariant, oldVariant);
                        }
                    }
                }
            }
        }
    }

    private void process(@NotNull final SageVariant mnv, @NotNull final SageVariant snv)
    {
        if(longerContainsShorter(snv, mnv))
        {
            snv.mixedGermlineImpact(mnv.mixedGermlineImpact());

            final Optional<GenomeRegion> maybeCodon = codon(snv.position());
            if(maybeCodon.filter(x -> keepMnv(x, snv.variant(), mnv.variant())).isPresent())
            {
                snv.filters().add(VariantVCF.DEDUP_FILTER);
            }
            else
            {
                mnv.filters().add(VariantVCF.DEDUP_FILTER);
            }
        }
    }

    @NotNull
    private Optional<GenomeRegion> codon(long position)
    {
        final Optional<HmfTranscriptRegion> maybeTranscript = mSelector.select(position);
        if(maybeTranscript.isPresent())
        {
            final List<GenomeRegion> codons = maybeTranscript.get().codonRangeAtGenomicPosition(position);
            for(GenomeRegion codon : codons)
            {
                if(position >= codon.start() && position <= codon.end())
                {
                    return Optional.of(codon);
                }
            }
        }
        return Optional.empty();
    }

    private static boolean isPassingSnv(@NotNull final SageVariant variant)
    {
        return variant.isPassing() && variant.isSnv();
    }

    private static boolean isMixedGermlineMnv(@NotNull final SageVariant variant)
    {
        return variant.isPassing() && variant.isMnv() && variant.mixedGermlineImpact() > 0;
    }

    static boolean keepMnv(@NotNull final GenomeRegion codon, @NotNull final VariantHotspot somaticSnv,
            @NotNull final VariantHotspot mixedMnv)
    {
        int snvCodonDifferences = codonDifferences(codon, somaticSnv);
        int mnvCodonDifferences = codonDifferences(codon, mixedMnv);

        return mnvCodonDifferences > snvCodonDifferences;
    }

    static int codonDifferences(@NotNull final GenomeRegion codon, @NotNull final VariantHotspot variant)
    {
        if(codon.start() > variant.end() || variant.position() > codon.end())
        {
            return 0;
        }

        int overlapStart = (int) Math.max(codon.start(), variant.position());
        int overlapEnd = (int) Math.min(codon.end(), variant.end());

        int difference = 0;
        for(int position = overlapStart; position <= overlapEnd; position++)
        {
            int variantIndex = (int) (position - variant.position());

            if(variant.ref().charAt(variantIndex) != variant.alt().charAt(variantIndex))
            {
                difference++;
            }
        }
        return difference;
    }
}
