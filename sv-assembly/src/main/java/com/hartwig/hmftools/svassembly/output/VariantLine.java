package com.hartwig.hmftools.svassembly.output;

import java.util.List;

import com.hartwig.hmftools.svassembly.models.AssemblyClassification;

import org.jetbrains.annotations.Nullable;

public class VariantLine
{
    @Nullable
    public final String LeftChromosome;
    public final int LeftPosition;
    @Nullable
    public final String RightChromosome;
    public final int RightPosition;
    @Nullable
    public final String LeftDescriptor;
    @Nullable
    public final String RightDescriptor;
    public final int MappingQuality;
    public final int GermlineSupport;
    public final int SomaticSupport;
    public final List<String> AssociatedAssemblies;
    public final AssemblyClassification Classification;
    public final List<String> SupportingReads;
    public final String Filters;

    public VariantLine(@Nullable final String leftChromosome, final int leftPosition, @Nullable final String rightChromosome,
            final int rightPosition, @Nullable final String leftDescriptor, @Nullable final String rightDescriptor,
            final int mappingQuality, final int germlineSupport, final int somaticSupport, final List<String> associatedAssemblies,
            final AssemblyClassification classification, final List<String> supportingReads, final String filters)
    {
        LeftChromosome = leftChromosome;
        LeftPosition = leftPosition;
        RightChromosome = rightChromosome;
        RightPosition = rightPosition;
        LeftDescriptor = leftDescriptor;
        RightDescriptor = rightDescriptor;
        MappingQuality = mappingQuality;
        GermlineSupport = germlineSupport;
        SomaticSupport = somaticSupport;
        AssociatedAssemblies = associatedAssemblies;
        Classification = classification;
        SupportingReads = supportingReads;
        Filters = filters;
    }
}
