package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// Delegate all method calls to the inner Variant
// This is created to allow VCF reading code to share a variant builder
public interface VariantDelegate extends Variant
{
    @NotNull
    Variant variant();

    @Override
    default int position()
    {
        return variant().position();
    }

    @Override
    default String chromosome()
    {
        return variant().chromosome();
    }

    @Override
    @NotNull
    default VariantType type()
	{
        return variant().type();
    }

    @Override
    @NotNull
    default String gene()
	{
        return variant().gene();
    }

    @Override
    @NotNull
    default String ref()
	{
        return variant().ref();
    }

    @Override
    @NotNull
    default String alt()
	{
        return variant().alt();
    }

    @Override
    @NotNull
    default AllelicDepth allelicDepth()
	{
        return variant().allelicDepth();
    }

    @Override
    @NotNull
    default String canonicalTranscript()
	{
        return variant().canonicalTranscript();
    }

    @Override
    @NotNull
    default String canonicalEffect()
	{
        return variant().canonicalEffect();
    }

    @Override
    @NotNull
    default CodingEffect canonicalCodingEffect()
	{
        return variant().canonicalCodingEffect();
    }

    @Override
    @NotNull
    default String canonicalHgvsCodingImpact()
	{
        return variant().canonicalHgvsCodingImpact();
    }

    @Override
    @NotNull
    default String canonicalHgvsProteinImpact()
	{
        return variant().canonicalHgvsProteinImpact();
    }

    @Override
    default double qual()
	{
        return variant().qual();
    }

    @Override
    default double mappability()
	{
        return variant().mappability();
    }

    @Override
    @NotNull
    default String filter()
	{
        return variant().filter();
    }

    @Override
    default int genesAffected()
	{
        return variant().genesAffected();
    }

    @Override
    default boolean spliceRegion()
	{
        return variant().spliceRegion();
    }

    @Override
    @NotNull
    default String otherReportedEffects()
	{
        return variant().otherReportedEffects();
    }

    @Override
    @NotNull
    default CodingEffect worstCodingEffect()
	{
        return variant().worstCodingEffect();
    }

    @Override
    @NotNull
    default Hotspot hotspot()
	{
        return variant().hotspot();
    }

    @Override
    default double adjustedCopyNumber()
	{
        return variant().adjustedCopyNumber();
    }

    @Override
    default double adjustedVAF()
	{
        return variant().adjustedVAF();
    }

    @Override
    default double minorAlleleCopyNumber()
	{
        return variant().minorAlleleCopyNumber();
    }

    @Override
    default double variantCopyNumber()
	{
        return variant().variantCopyNumber();
    }

    @Override
    default boolean biallelic()
	{
        return variant().biallelic();
    }

    @Override
    default boolean reported()
	{
        return variant().reported();
    }

    @Override
    @NotNull
    default GenotypeStatus genotypeStatus()
	{
        return variant().genotypeStatus();
    }

    @Override
    @NotNull
    default GermlineStatus germlineStatus()
	{
        return variant().germlineStatus();
    }

    @Override
    @NotNull
    default String trinucleotideContext()
	{
        return variant().trinucleotideContext();
    }

    @Override
    @NotNull
    default String microhomology()
	{
        return variant().microhomology();
    }

    @Override
    @NotNull
    default String repeatSequence()
	{
        return variant().repeatSequence();
    }

    @Override
    default int repeatCount()
	{
        return variant().repeatCount();
    }

    @Override
    @NotNull
    default VariantTier tier()
	{
        return variant().tier();
    }

    @Override
    @Nullable
    default AllelicDepth rnaDepth()
	{
        return variant().rnaDepth();
    }
    
    @Override
    @Nullable
    default List<String> reportableTranscripts()
	{
        return variant().reportableTranscripts();
    }

    @Override
    @Nullable
    default List<Integer> localPhaseSets()
    {
        return variant().localPhaseSets();
    }

    @Override
    @Nullable
    default String clinvarInfo()
	{
        return variant().clinvarInfo();
    }
}
