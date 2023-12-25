package com.hartwig.hmftools.esvee.processor;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.models.AlignedAssembly;
import com.hartwig.hmftools.esvee.models.Alignment;
import com.hartwig.hmftools.esvee.models.AssemblyClassification;
import com.hartwig.hmftools.esvee.models.AssemblyClassificationType;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.models.SupportedAssembly;
import com.hartwig.hmftools.esvee.util.ParallelMapper;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.util.SequenceUtil;

public class VariantCaller
{
    private final ExecutorService mExecutor;
    private final AnchorCigarFactory mAnchorFactory;

    public VariantCaller(final ExecutorService executor)
    {
        mExecutor = executor;
        mAnchorFactory = new AnchorCigarFactory();
    }

    public List<VariantCall> callVariants(final List<AlignedAssembly> assemblies)
    {
        return ParallelMapper.flatMap(mExecutor, assemblies, assembly -> callVariants(assembly));
    }

    public List<VariantCall> callVariants(final AlignedAssembly assembly)
    {
        final List<VariantCall> calls = new ArrayList<>();

        Alignment next = null;
        Alignment current = null;
        Alignment previous;
        for(Alignment alignment : assembly.getAlignmentBlocks())
        {
            previous = current;
            current = next;
            next = alignment;
            if(current != null && current.Chromosome.equals("-") && !calls.isEmpty())
            {
                return calls;
            }

            if(current == null || current.Chromosome.equals("-"))
                continue;
            @Nullable
            final Alignment nextValid = next.Chromosome.equals("-") ? null : next;

            @Nullable
            final var call = createCandidate(assembly, previous, current, nextValid);
            if(call != null)
                calls.add(call);
        }

        @Nullable
        final var last = createCandidate(assembly, current, next, null);
        if(last != null)
            calls.add(last);

        return calls;
    }

    @Nullable
    private VariantCall createCandidate(final AlignedAssembly assembly,
            @Nullable final Alignment previous, final Alignment current, @Nullable final Alignment next)
    {
        if(next == null && previous == null)
            return null; // Need at least 2 alignments to make a call

        @Nullable
        VariantCall result;
        if(next != null && previous != null)
        {
            if((result = tryCallTranslocationWithInsert(assembly, previous, current, next)) != null)
                return result;

            if((result = tryCallTranslocation(assembly, previous, current)) != null)
                return result;

            if((result = tryCallInsertDeleteOrDuplication(assembly, previous, current, next)) != null)
                return result;

            if((result = tryCallDeleteOrDuplication(assembly, previous, current)) != null)
                return result;

            return null;
        }
        if(previous == null)
        {
            if(next.isMapped() && current.isUnmapped() && current.Length > SvConstants.CALLERMINSIZETOCALL)
                return buildSingleEndedRight(assembly, current, next);
        }
        if(next == null)
        {
            if((result = tryCallTranslocation(assembly, previous, current)) != null)
                return result;
            if((result = tryCallDeleteOrDuplication(assembly, previous, current)) != null)
                return result;

            if(previous.isMapped() && current.isUnmapped() && current.Length > SvConstants.CALLERMINSIZETOCALL)
                return buildSingleEndedLeft(assembly, previous, current);
        }

        return null;
    }

    @Nullable
    private VariantCall tryCallTranslocationWithInsert(final AlignedAssembly assembly,
            final Alignment previous, final Alignment current, final Alignment next)
    {
        if(current.isMapped() || previous.isUnmapped() || next.isUnmapped())
            return null;
        if(next.Chromosome.equals(previous.Chromosome) && next.Inverted == previous.Inverted)
            return null;

        return buildTranslocation(assembly, previous, current, next);
    }

    @Nullable
    private VariantCall tryCallTranslocation(final AlignedAssembly assembly, final Alignment previous, final Alignment current)
    {
        if(current.isUnmapped() || previous.isUnmapped())
            return null;
        if(current.Chromosome.equals(previous.Chromosome) && current.Inverted == previous.Inverted)
            return null;

        return buildTranslocation(assembly, previous, null, current);
    }

    private String getAnchor(final AlignedAssembly assembly, final int position, final boolean invert)
    {
        final byte base = assembly.getBases()[position - 1];
        return String.valueOf((char) (invert ? SequenceUtil.complement(base) : base));
    }

    private VariantCall buildTranslocation(final AlignedAssembly assembly,
            final Alignment left, @Nullable final Alignment insert, final Alignment right)
    {
        final int leftEndSeqPosition = left.SequenceStartPosition + left.Length - 1;
        final int rightStartSeqPosition = right.SequenceStartPosition;
        final String leftAnchor = getAnchor(assembly, leftEndSeqPosition, false);
        final String rightAnchor = getAnchor(assembly, rightStartSeqPosition, false);
        final String insertSequence = insert == null
                ? ""
                : assembly.subsequence(insert.SequenceStartPosition - 1, insert.Length).getBasesString();
        final String leftInsertSequence = left.Inverted
                ? SequenceUtil.reverseComplement(leftAnchor + insertSequence)
                : leftAnchor + insertSequence;
        final String rightInsertSequence = right.Inverted
                ? SequenceUtil.reverseComplement(insertSequence + rightAnchor)
                : insertSequence + rightAnchor;

        final int leftEndPosition = endPosition(left, left.ReferenceStartPosition);
        final char leftDirection = left.Inverted ? '[' : ']';
        final String leftAsDestination = leftDirection + left.Chromosome + ":" + leftEndPosition + leftDirection;

        final int rightStartPosition = startPosition(right, right.ReferenceStartPosition);
        final char rightDirection = right.Inverted ? ']' : '[';
        final String rightAsDestination = rightDirection + right.Chromosome + ":" + rightStartPosition + rightDirection;

        final String leftDescriptor = left.Inverted
                ? rightAsDestination + leftInsertSequence
                : leftInsertSequence + rightAsDestination;
        final String rightDescriptor = right.Inverted
                ? rightInsertSequence + leftAsDestination
                : leftAsDestination + rightInsertSequence;

        final Support support = calculateSupport(assembly, left, right, insert);
        final AssemblyClassification classification = left.Chromosome.equals(right.Chromosome)
                ? new AssemblyClassification(AssemblyClassificationType.INVERSION, 0)
                : new AssemblyClassification(AssemblyClassificationType.TRANSLOCATION, 0);
        return build(left.Chromosome, leftEndPosition, right.Chromosome, rightStartPosition,
                leftDescriptor, rightDescriptor, assembly, left, right, insert, support, left.Quality, right.Quality,
                classification);
    }

    private int startPosition(final Alignment alignment, final int position)
    {
        return alignment.Inverted
                ? position + (alignment.Length - 1)
                : position;
    }

    private int endPosition(final Alignment alignment, final int position)
    {
        return alignment.Inverted
                ? position
                : position + (alignment.Length - 1);
    }

    @Nullable
    private VariantCall tryCallInsertDeleteOrDuplication(final AlignedAssembly assembly,
            final Alignment previous, @Nullable final Alignment current, final Alignment next)
    {
        if((current != null && current.isMapped()) || previous.isUnmapped() || next.isUnmapped())
            return null;
        if(!next.Chromosome.equals(previous.Chromosome) || next.Inverted != previous.Inverted)
            return null;

        // A delete happens when we move forward in the chromosome, a duplication happens if we move backwards
        final int previousEnd = previous.ReferenceStartPosition + previous.Length - 1;
        final int skippedBases = next.ReferenceStartPosition - previousEnd - 1;
        final int insertSize = current == null ? 0 : current.Length;
        if(insertSize + Math.abs(skippedBases) < SvConstants.CALLERMINSIZETOCALL)
            return null; // Not interested

        // We can either call it as an insert, a duplication, or a delete.
        if(current != null && insertSize >= Math.abs(skippedBases))
            return buildInsert(assembly, previous, current, next, insertSize);

        if(skippedBases > 0)
            return buildDelete(assembly, previous, current, next);
        else
            return buildDuplication(assembly, previous, current, next);
    }

    @Nullable
    private VariantCall tryCallDeleteOrDuplication(final AlignedAssembly assembly, final Alignment previous, final Alignment current)
    {
        return tryCallInsertDeleteOrDuplication(assembly, previous, null, current);
    }

    private VariantCall buildInsert(final AlignedAssembly assembly, final Alignment left, final Alignment insert, final Alignment right,
            final int insertSize)
    {
        final int leftEndSeqPosition = endPosition(left, left.SequenceStartPosition);
        final int rightStartSeqPosition = startPosition(right, right.SequenceStartPosition);
        final String leftAnchor = getAnchor(assembly, leftEndSeqPosition, left.Inverted);
        final String rightAnchor = getAnchor(assembly, rightStartSeqPosition, right.Inverted);
        final String insertSequence = assembly.subsequence(insert.SequenceStartPosition - 1, insert.Length).getBasesString();

        final int leftEndPosition = endPosition(left, left.ReferenceStartPosition);
        final int rightStartPosition = startPosition(right, right.ReferenceStartPosition);
        final String leftAsDestination = "]" + left.Chromosome + ":" + leftEndPosition + "]";
        final String rightAsDestination = "[" + right.Chromosome + ":" + rightStartPosition + "[";
        final String leftDescriptor = leftAnchor + insertSequence + rightAsDestination;
        final String rightDescriptor = leftAsDestination + insertSequence + rightAnchor;

        final Support support = calculateSupport(assembly, left, right, insert);
        return build(left.Chromosome, leftEndPosition, right.Chromosome, rightStartPosition,
                leftDescriptor, rightDescriptor, assembly, left, right, insert, support, left.Quality, right.Quality,
                new AssemblyClassification(AssemblyClassificationType.INSERT, insertSize));
    }

    private VariantCall buildDelete(final AlignedAssembly assembly,
            final Alignment left, @Nullable final Alignment insert, final Alignment right)
    {
        final int leftEndSeqPosition = endPosition(left, left.SequenceStartPosition);
        final int rightStartSeqPosition = startPosition(right, right.SequenceStartPosition);
        final String leftAnchor = getAnchor(assembly, leftEndSeqPosition, left.Inverted);
        final String rightAnchor = getAnchor(assembly, rightStartSeqPosition, right.Inverted);

        final int leftEndPosition = endPosition(left, left.ReferenceStartPosition);
        final int rightStartPosition = startPosition(right, right.ReferenceStartPosition);
        final String leftAsDestination = "]" + left.Chromosome + ":" + leftEndPosition + "]";
        final String rightAsDestination = "[" + right.Chromosome + ":" + rightStartPosition + "[";
        final String insertSequence = insert == null
                ? ""
                : assembly.subsequence(insert.SequenceStartPosition - 1, insert.Length).getBasesString();
        final String leftDescriptor = leftAnchor + insertSequence + rightAsDestination;
        final String rightDescriptor = leftAsDestination + insertSequence + rightAnchor;

        final Support support = calculateSupport(assembly, left, right, insert);
        return build(left.Chromosome, leftEndPosition, right.Chromosome, rightStartPosition,
                leftDescriptor, rightDescriptor, assembly, left, right, insert, support, left.Quality, right.Quality,
                new AssemblyClassification(AssemblyClassificationType.DELETION, rightStartPosition - leftEndPosition - 1));
    }

    private VariantCall buildDuplication(final AlignedAssembly assembly,
            final Alignment left, @Nullable final Alignment insert, final Alignment right)
    {
        final int leftEndSeqPosition = endPosition(left, left.SequenceStartPosition);
        final int rightStartSeqPosition = startPosition(right, right.SequenceStartPosition);
        final String leftAnchor = assembly.subsequence(leftEndSeqPosition - 1, 1).getBasesString();
        final String rightAnchor = assembly.subsequence(rightStartSeqPosition - 1, 1).getBasesString();

        final int leftEndPosition = endPosition(left, left.ReferenceStartPosition);
        final int rightStartPosition = startPosition(right, right.ReferenceStartPosition);
        final String leftAsDestination = "]" + left.Chromosome + ":" + leftEndPosition + "]";
        final String rightAsDestination = "[" + right.Chromosome + ":" + rightStartPosition + "[";
        final String insertSequence = insert == null
                ? ""
                : assembly.subsequence(insert.SequenceStartPosition - 1, insert.Length).getBasesString();
        final String leftDescriptor = leftAnchor + insertSequence + rightAsDestination;
        final String rightDescriptor = leftAsDestination + insertSequence + rightAnchor;

        final Support support = calculateSupport(assembly, left, right, insert);

        final int previousEndPosition = left.ReferenceStartPosition + left.Length;
        final int duplicationSize = previousEndPosition - rightStartPosition;
        final AssemblyClassification classification;
        if(duplicationSize < 6) // Can't call a duplication this short
            classification = new AssemblyClassification(AssemblyClassificationType.INSERT, duplicationSize);
        else
            classification = new AssemblyClassification(AssemblyClassificationType.DUPLICATION, duplicationSize);

        return build(left.Chromosome, leftEndPosition, right.Chromosome, rightStartPosition,
                leftDescriptor, rightDescriptor, assembly, left, right, insert, support, left.Quality, right.Quality, classification);
    }

    private VariantCall buildSingleEndedLeft(final AlignedAssembly assembly, final Alignment left, final Alignment insert)
    {
        final int leftEndSeqPosition = left.SequenceStartPosition + left.Length;
        final String leftAnchor = getAnchor(assembly, leftEndSeqPosition, left.Inverted);
        final String insertSequence = assembly.subsequence(leftEndSeqPosition - 1, assembly.getLength() - leftEndSeqPosition - 1).getBasesString();
        final String leftDescriptor = leftAnchor + insertSequence + ".";

        final Support support = calculateSupport(assembly, left, null, insert);
        return build(left.Chromosome, endPosition(left, left.ReferenceStartPosition), null, 0,
                leftDescriptor, null, assembly, left, null, insert, support, left.Quality, 0,
                new AssemblyClassification(AssemblyClassificationType.UNKNOWN, 0));
    }

    private VariantCall buildSingleEndedRight(final AlignedAssembly assembly, final Alignment insert, final Alignment right)
    {
        final String rightAnchor = getAnchor(assembly, right.SequenceStartPosition - 1, right.Inverted);
        final String insertSequence = assembly.subsequence(0, right.SequenceStartPosition - 1).getBasesString();
        final String rightDescriptor = "." + insertSequence + rightAnchor;

        final Support support = calculateSupport(assembly, null, right, insert);
        return build(null, 0,
                right.Chromosome, startPosition(right, right.ReferenceStartPosition),
                null,  rightDescriptor, assembly, null, right, insert, support,
                0, right.Quality, new AssemblyClassification(AssemblyClassificationType.UNKNOWN, 0));
    }

    private VariantCall build(
            @Nullable final String leftChromosome, final int leftPosition,
            @Nullable final String rightChromosome, final int rightPosition,
            @Nullable final String leftDescriptor, @Nullable final String rightDescriptor,
            final AlignedAssembly assembly, @Nullable final Alignment left, @Nullable final Alignment right,
            @Nullable final Alignment insert, final Support support,
            final int leftQuality, final int rightQuality, final AssemblyClassification classification)
    {
        if (leftChromosome != null && rightChromosome != null)
        {
            if (leftQuality < 30 && rightQuality >= 30)
                return buildSingleEndedRight(assembly, Objects.requireNonNullElse(insert, left), Objects.requireNonNull(right));
            else if (leftQuality >= 30 && rightQuality < 30)
                return buildSingleEndedLeft(assembly, Objects.requireNonNull(left), Objects.requireNonNullElse(insert, right));
        }

        final List<VariantCall.SampleSupport> sampleSupport = partitionSupport(support);
        final var anchors = mAnchorFactory.anchorCigar(assembly, left, right);
        final int leftOffset = left != null ? left.SequenceStartPosition + left.Length - 1 : 0;
        final int leftOverhang = left != null ? calculateLeftOverhang(assembly, leftOffset, support) : 0;
        final int rightOffset = right != null ? right.SequenceStartPosition : 0;
        final int rightOverhang = right != null ? calculateRightOverhang(assembly, rightOffset, support) : 0;
        final VariantCall.VariantAssembly variantAssembly = VariantCall.VariantAssembly.create(assembly,
                anchors.getLeft(), leftOffset, leftOverhang,
                anchors.getRight(), rightOffset, rightOverhang);
        return VariantCall.create(leftChromosome, leftPosition, rightChromosome, rightPosition, leftDescriptor, rightDescriptor,
                Set.of(1), Set.of(variantAssembly), leftQuality, rightQuality, sampleSupport, classification);
    }

    private int calculateLeftOverhang(final SupportedAssembly assembly, final int leftOffset, final Support support)
    {
        int maxOverhang = 0;
        for(Record record : support.SplitReads)
        {
            final int supportStartOffset = assembly.getSupportIndex(record) + 1;
            final int supportEndOffset = assembly.getSupportIndex(record) + record.getLength();
            if (supportStartOffset > leftOffset || supportEndOffset < leftOffset)
                continue; // Doesn't cross left
            final int overhang = leftOffset - supportStartOffset + 1;
            maxOverhang = Math.max(overhang, maxOverhang);
        }
        return maxOverhang;
    }

    private int calculateRightOverhang(final SupportedAssembly assembly, final int rightOffset, final Support support)
    {
        int maxOverhang = 0;
        for(Record record : support.SplitReads)
        {
            final int supportStartOffset = assembly.getSupportIndex(record) + 1;
            final int supportEndOffset = assembly.getSupportIndex(record) + record.getLength();
            if (supportEndOffset < rightOffset || supportStartOffset > rightOffset)
                continue; // Doesn't cross right
            final int overhang = supportEndOffset - rightOffset + 1;
            maxOverhang = Math.max(overhang, maxOverhang);
        }
        return maxOverhang;
    }

    private Support calculateSupport(final AlignedAssembly assembly, @Nullable final Alignment left, @Nullable final Alignment right,
            @Nullable final Alignment insert)
    {
        final var splitReads = findSplitReads(assembly, left, right, insert);
        final var discordantPairs = findDiscordantPairs(assembly, left, right, insert, splitReads);
        return new Support(splitReads, discordantPairs);
    }

    private Set<Record> findSplitReads(final AlignedAssembly assembly, @Nullable final Alignment left, @Nullable final Alignment right,
            @Nullable final Alignment insert)
    {
        final Set<Record> splitReads = new HashSet<>();
        for(Map.Entry<Record, Integer> entry : assembly.getSupport())
        {
            final Record record = entry.getKey();
            final int supportLeft = entry.getValue();
            final int supportRight = supportLeft + record.getLength();

            if(left != null)
            {
                final int leftAlignmentRight = left.SequenceStartPosition + left.Length;
                if(supportLeft <= leftAlignmentRight && supportRight > leftAlignmentRight)
                    splitReads.add(record);
            }

            if(right != null)
            {
                final int rightAlignmentLeft = right.SequenceStartPosition;
                if(supportLeft < rightAlignmentLeft && supportRight >= rightAlignmentLeft)
                    splitReads.add(record);
            }

            if(insert != null)
            {
                final int alignmentLeft = insert.SequenceStartPosition;
                final int alignmentRight = insert.SequenceStartPosition + insert.Length;
                if(supportLeft <= alignmentRight && supportRight >= alignmentLeft)
                    splitReads.add(record);
            }
        }

        // Add mates of anything split, if supporting
        final Set<String> splitReadFragments = splitReads.stream()
                .map(Record::getName)
                .collect(Collectors.toSet());
        assembly.getSupportRecords().stream()
                .filter(record -> splitReadFragments.contains(record.getName()))
                .forEach(splitReads::add);

        return splitReads;
    }

    private Set<Record> findDiscordantPairs(final AlignedAssembly assembly, @Nullable final Alignment left, @Nullable final Alignment right,
            @Nullable final Alignment insert, final Set<Record> splitReads)
    {
        if (left == null || right == null)
            return new HashSet<>();

        final Set<Record> support = intersectionSupport(
                supportingReadsLeft(assembly, left),
                supportingReadsRight(assembly, right),
                supportingReads(assembly, insert));

        final Set<String> splitReadFragments = splitReads.stream()
                .map(Record::getName)
                .collect(Collectors.toSet());
        support.removeIf(r -> splitReadFragments.contains(r.getName()));

        final Set<String> discordantFragments = support.stream()
                .filter(record -> record.isDiscordant(SvConstants.DISCORDANTPAIRFRAGMENTLENGTH)
                        || record.isUnmapped() || record.isMateUnmapped())
                .map(Record::getName)
                .collect(Collectors.toSet());
        support.removeIf(r -> !discordantFragments.contains(r.getName()));

        return support;
    }

    private Set<Record> intersectionSupport(final Set<Record> leftSupport, final Set<Record> rightSupport, final Set<Record> insertSupport)
    {
        final Set<Record> support = new HashSet<>();
        final var insertSupportFragments = insertSupport.stream().map(Record::getName).collect(Collectors.toSet());
        support.addAll(insertSupport);

        final var leftSupportFragments = leftSupport.stream().map(Record::getName).collect(Collectors.toSet());
        leftSupportFragments.addAll(insertSupportFragments);

        final var rightSupportFragments = rightSupport.stream().map(Record::getName).collect(Collectors.toSet());
        rightSupportFragments.addAll(insertSupportFragments);

        final Set<String> bothSideFragments = new HashSet<>(leftSupportFragments);
        bothSideFragments.retainAll(rightSupportFragments);

        leftSupport.stream().filter(r -> bothSideFragments.contains(r.getName())).forEach(support::add);
        rightSupport.stream().filter(r -> bothSideFragments.contains(r.getName())).forEach(support::add);

        return support;
    }

    private Set<Record> supportingReads(final SupportedAssembly assembly, @Nullable final Alignment alignment)
    {
        if(alignment == null)
            return Set.of();

        final Set<Record> support = new HashSet<>();
        for(Map.Entry<Record, Integer> entry : assembly.getSupport())
        {
            final Record record = entry.getKey();
            final int supportLeft = entry.getValue();
            final int supportRight = supportLeft + record.getLength();

            final int alignmentLeft = alignment.SequenceStartPosition;
            final int alignmentRight = alignment.SequenceStartPosition + alignment.Length;
            if(supportLeft <= alignmentRight && supportRight >= alignmentLeft)
                support.add(record);
        }

        return support;
    }

    private Set<Record> supportingReadsLeft(final AlignedAssembly assembly, @Nullable final Alignment stopAt)
    {
        if(stopAt == null)
            return Set.of();

        final Set<Record> support = new HashSet<>();
        for(Alignment alignment : assembly.getAlignmentBlocks())
        {
            support.addAll(supportingReads(assembly, alignment));
            if(alignment == stopAt)
                break;
        }
        return support;
    }

    private Set<Record> supportingReadsRight(final AlignedAssembly assembly, @Nullable final Alignment startAt)
    {
        if(startAt == null)
            return Set.of();

        final Set<Record> support = new HashSet<>();

        int alignmentIndex = 0;
        for(; alignmentIndex < assembly.getAlignmentBlocks().size(); alignmentIndex++)
            if(assembly.getAlignmentBlocks().get(alignmentIndex) == startAt)
                break;
        for(; alignmentIndex < assembly.getAlignmentBlocks().size(); alignmentIndex++)
            support.addAll(supportingReads(assembly, assembly.getAlignmentBlocks().get(alignmentIndex)));

        return support;
    }

    private List<VariantCall.SampleSupport> partitionSupport(final Support support)
    {
        final Map<String, Pair<Set<Record>, Set<Record>>> bySample = new LinkedHashMap<>();

        Stream.concat(support.DiscordantSupport.stream(), support.SplitReads.stream())
                .map(Record::sampleName)
                .distinct()
                .sorted()
                .forEach(sampleName -> bySample.put(sampleName, Pair.of(new HashSet<>(), new HashSet<>())));

        for (final Record record : support.SplitReads)
            bySample.get(record.sampleName()).getLeft().add(record);
        for (final Record record : support.DiscordantSupport)
            bySample.get(record.sampleName()).getRight().add(record);

        final List<VariantCall.SampleSupport> sampleSupport = new ArrayList<>();
        for(Map.Entry<String, Pair<Set<Record>, Set<Record>>> entry : bySample.entrySet())
        {
            final String sampleName = entry.getKey();
            final Set<Record> splitReads = entry.getValue().getLeft();
            final Set<Record> discordantReads = entry.getValue().getRight();

            final boolean isGermline = Stream.concat(splitReads.stream(), discordantReads.stream()).anyMatch(Record::isGermline);
            final int quality = splitReads.size() + discordantReads.size();

            sampleSupport.add(new VariantCall.SampleSupport(sampleName, isGermline, quality, splitReads, discordantReads));
        }

        return sampleSupport;
    }

    private static class Support
    {
        public final Set<Record> SplitReads;
        public final Set<Record> DiscordantSupport;

        private Support(final Set<Record> splitReads, final Set<Record> discordantSupport)
        {
            SplitReads = splitReads;
            DiscordantSupport = discordantSupport;
        }
    }
}
