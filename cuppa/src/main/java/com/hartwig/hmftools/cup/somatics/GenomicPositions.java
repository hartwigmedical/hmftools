package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.sigs.SnvSigUtils.variantContext;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.AID_APOBEC_TRINUCLEOTIDE_CONTEXTS;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.ALL;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.FALSE_ONLY;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.TRUE_ONLY;

import java.util.List;

import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.variant.VariantType;

public final class GenomicPositions
{
    public static void extractPositionFrequencyCounts(
            final List<SomaticVariant> variants, final PositionFrequencies posFrequencies, AidApobecStatus aidApobecStatus)
    {
        posFrequencies.clear();

        for(final SomaticVariant variant : variants)
        {
            if(variant.Type != VariantType.SNP)
                continue;

            // exclude male chromosome since is then unhelpful for multi-gender cancer types
            if(variant.Chromosome.equals("Y") || variant.Chromosome.equals("chrY"))
                continue;

            if(variant.TrinucleotideContext.contains("N"))
                continue;

            if(aidApobecStatus != ALL)
            {
                String bucketName = variantContext(variant.Ref, variant.Alt, variant.TrinucleotideContext);
                boolean isAA = AID_APOBEC_TRINUCLEOTIDE_CONTEXTS.contains(bucketName);

                if((aidApobecStatus == TRUE_ONLY && !isAA) || (aidApobecStatus == FALSE_ONLY && isAA))
                    continue;
            }

            if(!posFrequencies.isValidChromosome(variant.Chromosome))
            {
                CUP_LOGGER.warn("variant chr({}) position({}) cannot map to genomic position", variant.Chromosome, variant.Position);
                continue;
            }

            posFrequencies.addPosition(variant.Chromosome, variant.Position);
        }
    }
}
