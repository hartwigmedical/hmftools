package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.variant.PurpleVcfTags;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.AllelicDepth;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class GermlinePurityEnrichment
{
    private final String mVersion;
    private final String mTumorSample;
    private final String mReferenceSample;
    private final PurityAdjuster mPurityAdjuster;
    private final GenomeRegionSelector<PurpleCopyNumber> mCopyNumberSelector;

    public GermlinePurityEnrichment(
            final String purpleVersion, final String tumorSample, final String referenceSample,
            @Nullable  final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers)
    {
        mVersion = purpleVersion;
        mTumorSample = tumorSample;
        mReferenceSample = referenceSample;
        mPurityAdjuster = purityAdjuster;
        mCopyNumberSelector = GenomeRegionSelectorFactory.createImproved(Multimaps.fromRegions(copyNumbers));
    }

    public void processVariant(final VariantContext variant)
    {
        if(mPurityAdjuster == null)
        {
            variant.getCommonInfo().putAttribute(PURPLE_VARIANT_CN, 1);
            variant.getCommonInfo().putAttribute(PURPLE_CN, 2);
            variant.getCommonInfo().putAttribute(PURPLE_AF, 0.5);
            variant.getCommonInfo().putAttribute(PURPLE_MINOR_ALLELE_CN_INFO, 1);
            variant.getCommonInfo().putAttribute(PURPLE_BIALLELIC_FLAG, false);
            return;
        }

        final Genotype tumorGenotype = variant.getGenotype(mTumorSample);
        final Genotype normalGenotype = variant.getGenotype(mReferenceSample);
        if(tumorGenotype != null && normalGenotype != null && tumorGenotype.hasAD() && HumanChromosome.contains(variant.getContig()))
        {
            final GenomePosition position = GenomePositions.create(variant.getContig(), variant.getStart());
            final AllelicDepth tumorDepth = AllelicDepth.fromGenotype(tumorGenotype);
            final GenotypeStatus germlineGenotype = GenotypeStatus.fromGenotype(normalGenotype);

            Optional<PurpleCopyNumber> optionalPurpleCopyNumber = mCopyNumberSelector.select(position);
            if(optionalPurpleCopyNumber.isPresent())
            {
                PurpleCopyNumber purpleCopyNumber = optionalPurpleCopyNumber.get();
                double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
                double vaf = germlineGenotype == GenotypeStatus.HOM_ALT ? 1.0 : vaf(germlineGenotype, purpleCopyNumber, tumorDepth);
                double variantCopyNumber = Math.max(0, vaf * copyNumber);
                boolean biallelic = Doubles.lessOrEqual(copyNumber, 0) || Doubles.greaterOrEqual(variantCopyNumber, copyNumber - 0.5);

                variant.getCommonInfo().putAttribute(PURPLE_VARIANT_CN, variantCopyNumber);
                variant.getCommonInfo().putAttribute(PURPLE_CN, copyNumber);
                variant.getCommonInfo().putAttribute(PURPLE_AF, String.format("%.4f", vaf));
                variant.getCommonInfo().putAttribute(PURPLE_MINOR_ALLELE_CN_INFO, purpleCopyNumber.minorAlleleCopyNumber());
                variant.getCommonInfo().putAttribute(PURPLE_BIALLELIC_FLAG, biallelic);
            }
        }
    }

    private double vaf(final GenotypeStatus germlineGenotype, PurpleCopyNumber purpleCopyNumber, final AllelicDepth tumorDepth)
    {
        if(tumorDepth.TotalReadCount == 0 || tumorDepth.AlleleReadCount == 0)
            return 0;

        if(Doubles.lessOrEqual(purpleCopyNumber.averageTumorCopyNumber(), 0.001))
            return 0;

        double rawAF = tumorDepth.alleleFrequency();
        double constrainedCopyNumber = Math.max(0.001, purpleCopyNumber.averageTumorCopyNumber());

        if(germlineGenotype.equals(GenotypeStatus.HET))
        {
            return mPurityAdjuster.purityAdjustedVAFWithHeterozygousNormal(purpleCopyNumber.chromosome(), constrainedCopyNumber, rawAF);
        }
        else
        {
            return mPurityAdjuster.purityAdjustedVAFWithHomozygousNormal(purpleCopyNumber.chromosome(), constrainedCopyNumber, rawAF);
        }
    }

    public void enrichHeader(final VCFHeader header)
    {
        PurpleVcfTags.addGermlineHeader(mVersion, header);
    }
}
