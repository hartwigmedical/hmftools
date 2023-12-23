package com.hartwig.hmftools.esvee.processor;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.models.AlignedAssembly;
import com.hartwig.hmftools.esvee.models.AssemblyClassification;
import com.hartwig.hmftools.esvee.models.AssemblyClassificationType;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.util.NaturalSortComparator;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;

public class VariantCall
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
            Classification = new AssemblyClassification(AssemblyClassificationType.UNKNOWN, 0);
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

    public String leftRef(final Context context)
    {
        if(LeftChromosome == null)
            return "N";
        return context.ReferenceGenome.getBaseString(LeftChromosome, LeftPosition, LeftPosition);
    }

    public String rightRef(final Context context)
    {
        if(RightChromosome == null)
            return "N";
        return context.ReferenceGenome.getBaseString(RightChromosome, RightPosition, RightPosition);
    }

    public int quality()
    {
        final int mappingQuality = LeftChromosome != null && RightChromosome != null
                ? Math.min(LeftMappingQuality, RightMappingQuality)
                : LeftChromosome == null ? RightMappingQuality : LeftMappingQuality;
        final int supportQuality = (int) (Math.log(germlineSupport() + somaticSupport()) * 2);
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

    public int germlineSupport()
    {
        return mSampleSupport.stream()
                .filter(SampleSupport::isGermline)
                .mapToInt(SampleSupport::totalSupportFragmentCount)
                .sum();
    }

    public int somaticSupport()
    {
        return mSampleSupport.stream()
                .filter(sampleSupport -> !sampleSupport.isGermline())
                .mapToInt(SampleSupport::totalSupportFragmentCount)
                .sum();
    }

    public Set<String> supportingFragments()
    {
        return mSampleSupport.stream()
                .flatMap(support -> Stream.concat(support.splitReads().stream(), support.discordantReads().stream()))
                .map(Record::getName)
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

    public static class SampleSupport
    {
        private final String mSampleName;
        private final boolean mIsGermline;
        private final int mQuality;
        private final Set<Record> mSplitReads;
        private final Set<Record> mDiscordantReads;
        private final int mSplitReadFragmentCount;
        private final int mDiscordantPairFragmentCount;

        public SampleSupport(final String sampleName, final boolean isGermline,
                final int quality, final Set<Record> splitReads, final Set<Record> discordantReads)
        {
            mSampleName = sampleName;
            mIsGermline = isGermline;
            mQuality = quality;
            mSplitReads = splitReads;
            mDiscordantReads = discordantReads;

            mSplitReadFragmentCount = (int) splitReads.stream().map(Record::getName).distinct().count();
            mDiscordantPairFragmentCount = (int) discordantReads.stream().map(Record::getName).distinct().count();
        }

        public String sampleName()
        {
            return mSampleName;
        }

        public boolean isGermline()
        {
            return mIsGermline;
        }

        public int quality()
        {
            return mQuality;
        }

        public Set<Record> splitReads()
        {
            return mSplitReads;
        }

        public Set<Record> discordantReads()
        {
            return mDiscordantReads;
        }

        public int totalSupportFragmentCount()
        {
            return mSplitReadFragmentCount + mDiscordantPairFragmentCount;
        }

        public int splitReadFragmentCount()
        {
            return mSplitReadFragmentCount;
        }

        public int discordantPairFragmentCount()
        {
            return mDiscordantPairFragmentCount;
        }
    }

    public static class VariantAssembly
    {
        public final AlignedAssembly Assembly;
        @Nullable
        public final Cigar LeftAnchorCigar;
        @Nullable
        public final Cigar RightAnchorCigar;
        public final int LeftCigarLength;
        public final int RightCigarLength;
        /**
         * Position in assembly
         */
        public final int LeftPosition;
        /**
         * Position in assembly
         */
        public final int RightPosition;
        public final int LeftOverhang;
        public final int RightOverhang;

        public VariantAssembly(final AlignedAssembly assembly,
                @Nullable final Cigar leftAnchorCigar,
                final int leftCigarLength,
                final int leftPosition,
                final int leftOverhang,
                @Nullable final Cigar rightAnchorCigar,
                final int rightCigarLength,
                final int rightPosition,
                final int rightOverhang)
        {
            Assembly = assembly;
            LeftAnchorCigar = leftAnchorCigar;
            LeftCigarLength = leftCigarLength;
            LeftPosition = leftPosition;
            LeftOverhang = leftOverhang;
            RightAnchorCigar = rightAnchorCigar;
            RightCigarLength = rightCigarLength;
            RightPosition = rightPosition;
            RightOverhang = rightOverhang;
        }

        public static VariantAssembly create(final AlignedAssembly assembly,
                @Nullable final Pair<Cigar, Integer> leftAnchor,
                final int leftPosition, final int leftOverhang,
                @Nullable final Pair<Cigar, Integer> rightAnchor,
                final int rightPosition, final int rightOverhang)
        {
            return new VariantAssembly(
                    assembly,
                    leftAnchor == null ? null : leftAnchor.getKey(),
                    leftAnchor == null ? 0 : leftAnchor.getValue(),
                    leftPosition, leftOverhang,
                    rightAnchor == null ? null : rightAnchor.getKey(),
                    rightAnchor == null ? 0 : rightAnchor.getValue(),
                    rightPosition, rightOverhang
            );
        }

        public VariantAssembly reverse()
        {
            return new VariantAssembly(Assembly,
                    RightAnchorCigar,
                    RightCigarLength,
                    RightPosition,
                    RightOverhang,
                    LeftAnchorCigar,
                    LeftCigarLength,
                    LeftPosition,
                    LeftOverhang);
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
                return true;
            else if(o == null || getClass() != o.getClass())
                return false;

            final VariantAssembly that = (VariantAssembly) o;
            return LeftCigarLength == that.LeftCigarLength && LeftPosition == that.LeftPosition
                    && RightCigarLength == that.RightCigarLength && RightPosition == that.RightPosition
                    && Assembly.equals(that.Assembly)
                    && Objects.equals(LeftAnchorCigar, that.LeftAnchorCigar) && Objects.equals(RightAnchorCigar, that.RightAnchorCigar);
        }

        @Override
        public int hashCode()
        {
            int result = Assembly.hashCode();
            result = 31 * result + (LeftAnchorCigar != null ? LeftAnchorCigar.hashCode() : 0);
            result = 31 * result + (RightAnchorCigar != null ? RightAnchorCigar.hashCode() : 0);
            result = 31 * result + LeftCigarLength;
            result = 31 * result + RightCigarLength;
            result = 31 * result + LeftPosition;
            result = 31 * result + RightPosition;
            return result;
        }
    }
}
