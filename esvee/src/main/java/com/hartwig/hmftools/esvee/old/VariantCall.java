package com.hartwig.hmftools.esvee.old;

import java.util.LinkedHashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.read.Read;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

// refactored class will just be called Variant
public class VariantCall
{
    /*
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

    public final Set<Integer> PhaseSets;
    public final int LeftMappingQuality;
    public final int RightMappingQuality;

    public final AssemblyClassification Classification;

    private final Set<VariantAssembly> mAssociatedAssemblies;
    private final List<SampleSupport> mSampleSupport;

    public static VariantCall create(
            @Nullable final String leftChromosome, final int leftPosition,
            @Nullable final String rightChromosome, final int rightPosition,
            @Nullable final String leftDescriptor, @Nullable final String rightDescriptor, final Set<Integer> phaseSets,
            final Set<VariantAssembly> associatedAssemblies, final int leftMappingQuality, final int rightMappingQuality,
            final List<SampleSupport> sampleSupport, final AssemblyClassification classification)
    {
        if(leftChromosome != null && rightChromosome != null)
        {
            int compare = NaturalSortComparator.INSTANCE.compare(leftChromosome, rightChromosome);
            if(compare == 0)
                compare = Integer.compare(leftPosition, rightPosition);
            if(compare > 0)
            {
                final Set<VariantAssembly> swappedAssemblies = associatedAssemblies.stream()
                        .map(VariantAssembly::reverse)
                        .collect(Collectors.toCollection(LinkedHashSet::new));
                return create(rightChromosome, rightPosition,
                        leftChromosome, leftPosition,
                        rightDescriptor, leftDescriptor,
                        phaseSets, swappedAssemblies, rightMappingQuality, leftMappingQuality,
                        sampleSupport, classification);
            }
        }
        if(leftChromosome == null && rightChromosome != null)
        {
            final Set<VariantAssembly> swappedAssemblies = associatedAssemblies.stream()
                    .map(VariantAssembly::reverse)
                    .collect(Collectors.toCollection(LinkedHashSet::new));
            return create(rightChromosome, rightPosition,
                    null, 0,
                    rightDescriptor, null,
                    phaseSets, swappedAssemblies, rightMappingQuality, leftMappingQuality,
                    sampleSupport, classification);
        }
        if(leftChromosome != null && leftMappingQuality < 30 && rightMappingQuality > 30)
        {
            // Remove left
            return create(null, 0,
                    rightChromosome, rightPosition,
                    null, rightDescriptor,
                    phaseSets, associatedAssemblies, 0, rightMappingQuality,
                    sampleSupport, classification);
        }
        else if(rightChromosome != null && rightMappingQuality < 30 && leftMappingQuality > 30)
        {
            // Remove right
            return create(leftChromosome, leftPosition,
                    null, 0,
                    leftDescriptor, null,
                    phaseSets, associatedAssemblies, leftMappingQuality, 0,
                    sampleSupport, classification);
        }

        return new VariantCall(leftChromosome, leftPosition, rightChromosome, rightPosition,
                leftDescriptor, rightDescriptor, phaseSets, associatedAssemblies,
                leftMappingQuality, rightMappingQuality, sampleSupport, classification);
    }

    protected VariantCall(
            @Nullable final String leftChromosome, final int leftPosition,
            @Nullable final String rightChromosome, final int rightPosition,
            @Nullable final String leftDescriptor, @Nullable final String rightDescriptor, final Set<Integer> phaseSets,
            final Set<VariantAssembly> associatedAssemblies, final int leftMappingQuality, final int rightMappingQuality,
            final List<SampleSupport> sampleSupport, final AssemblyClassification classification)
    {
        assert leftChromosome != null || rightChromosome != null;
        LeftChromosome = leftChromosome;
        LeftPosition = leftPosition;
        LeftDescriptor = leftDescriptor;
        LeftMappingQuality = leftMappingQuality;
        RightChromosome = rightChromosome;
        RightPosition = rightPosition;
        RightDescriptor = rightDescriptor;
        RightMappingQuality = rightMappingQuality;

        if(LeftChromosome == null || RightChromosome == null)
            Classification = new AssemblyClassification(StructuralVariantType.SGL, 0);
        else
            Classification = classification;

        PhaseSets = phaseSets;
        mAssociatedAssemblies = associatedAssemblies;
        mSampleSupport = sampleSupport;
    }

    public boolean isGermline()
    {
        return germlineSupport() > 0;
    }

    public boolean isSingleSided()
    {
        return LeftChromosome == null || RightChromosome == null;
    }

    public String leftRef(final RefGenomeInterface refGenome)
    {
        // CHECK: use of N for SGLs?
        if(LeftChromosome == null)
            return "N";

        return refGenome.getBaseString(LeftChromosome, LeftPosition, LeftPosition);
    }

    public String rightRef(final RefGenomeInterface refGenome)
    {
        if(RightChromosome == null)
            return "N";

        return refGenome.getBaseString(RightChromosome, RightPosition, RightPosition);
    }

    public int quality()
    {
        int mappingQuality = LeftChromosome != null && RightChromosome != null
                ? Math.min(LeftMappingQuality, RightMappingQuality)
                : LeftChromosome == null ? RightMappingQuality : LeftMappingQuality;
        // final int supportQuality = (int) (Math.log(germlineSupport() + somaticSupport()) * 2);
        int supportQuality = (int) (Math.log(germlineSupport() + somaticSupport()) * 2);
        return mappingQuality + supportQuality;
    }



    @SuppressWarnings("DataFlowIssue")
    public String name()
    {
        if(isSingleSided())
        {
            final String baseName = LeftChromosome == null
                    ? RightChromosome + "_" + RightPosition
                    : LeftChromosome + "_" + LeftPosition;
            final String hash = String.valueOf(LeftDescriptor == null
                    ? RightDescriptor.hashCode()
                    : LeftDescriptor.hashCode());
            return baseName + "_" + hash;
        }
        final String leftName = LeftChromosome + "_" + LeftPosition;
        final String rightName = RightChromosome + "_" + RightPosition;
        final String hash = String.valueOf(LeftDescriptor.hashCode() ^ RightDescriptor.hashCode());
        return leftName + "_" + rightName + "_" + hash;
    }

    @SuppressWarnings("DataFlowIssue")
    public String compactName()
    {
        if(isSingleSided())
        {
            final String baseName = LeftChromosome == null
                    ? RightChromosome + Integer.toString(RightPosition, 36)
                    : LeftChromosome + Integer.toString(LeftPosition, 36);
            final String hash = Integer.toString(LeftDescriptor == null
                    ? RightDescriptor.hashCode()
                    : LeftDescriptor.hashCode(), 36);
            return baseName + hash;
        }
        final String leftName = LeftChromosome + Integer.toString(LeftPosition, 36);
        final String rightName = RightChromosome + Integer.toString(RightPosition, 36);
        final String hash = Integer.toString((LeftDescriptor.hashCode() * 13) ^ RightDescriptor.hashCode(), 36);
        return leftName + rightName + hash;
    }

    public int supportCount()
    {
        return mSampleSupport.stream()
                .mapToInt(SampleSupport::totalSupportFragmentCount)
                .sum();
    }

    public int germlineSupport()
    {
        return supportCount();
    }

    public int somaticSupport()
    {
        return supportCount();
    }

    public Set<String> supportingFragments()
    {
        return mSampleSupport.stream()
                .flatMap(support -> Stream.concat(support.splitReads().stream(), support.discordantReads().stream()))
                .map(Read::getName)
                .collect(Collectors.toSet());
    }

    public List<SampleSupport> sampleSupport()
    {
        return mSampleSupport;
    }

    public int splitReadSupport()
    {
        return sampleSupport().stream()
                .mapToInt(SampleSupport::splitReadFragmentCount)
                .sum();
    }

    public int discordantSupport()
    {
        return sampleSupport().stream()
                .mapToInt(SampleSupport::discordantPairFragmentCount)
                .sum();
    }

    public Set<AlignedAssembly> associatedAssemblies()
    {
        return mAssociatedAssemblies.stream()
                .map(variantAssembly -> variantAssembly.Assembly)
                .collect(Collectors.toCollection(LinkedHashSet::new));
    }

    public Set<VariantAssembly> variantAssemblies()
    {
        return mAssociatedAssemblies;
    }

    public int overhang()
    {
        if(isSingleSided())
        {
            return variantAssemblies().stream()
                    .mapToInt(assembly -> assembly.LeftOverhang)
                    .max()
                    .orElse(0);
        }
        else
        {
            final Pair<Integer, Integer> overhang = variantAssemblies().stream()
                    .map(assembly -> Pair.of(assembly.LeftOverhang, assembly.RightOverhang))
                    .reduce((l, r) -> Pair.of(Math.max(l.getLeft(), r.getLeft()), Math.max(l.getRight(), r.getRight())))
                    .orElse(Pair.of(0, 0));
            return Math.min(overhang.getLeft(), overhang.getRight());
        }
    }

    @Override
    public String toString()
    {
        return "{\n\t\"Left\": \"" + LeftChromosome + "@" + LeftPosition + "\",\n\t"
                + "\"LeftDescriptor\": \"" + LeftDescriptor + "\",\n\t"
                + "\"Right\": \"" + RightChromosome + "@" + RightPosition + "\",\n\t"
                + "\"RightDescriptor\": \"" + RightDescriptor + "\",\n\t"
                + "\"Support\": " + supportingFragments().size() + "\n}";
    }
    */
}
