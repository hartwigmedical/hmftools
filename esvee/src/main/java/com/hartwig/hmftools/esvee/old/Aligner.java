package com.hartwig.hmftools.esvee.old;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.TASK_LOG_COUNT;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.CigarUtils;
import com.hartwig.hmftools.esvee.variant.HomologySlider;
import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.common.ThreadTask;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFlag;

// @ThreadSafe
public class Aligner extends ThreadTask implements AutoCloseable
{
    private final Queue<GappedAssembly> mSecondaryAssemblyQueue;
    private final int mSecondaryAssemblyCount;

    private final SvConfig mConfig;
    private final BwaMemIndex mIndex;

    // CHECK: what benefit if any is there is using a ThreadLocal here?
    private final ThreadLocal<BwaMemAligner> mAligners;

    private final HomologySlider mHomologySlider;

    private final List<AlignedAssembly> mAlignedAssemblies;

    public static List<Aligner> createThreadTasks(
            final List<GappedAssembly> secondaryAssemblies, final SvConfig config, final int taskCount, final List<Thread> threadTasks)
    {
        List<Aligner> alignerTasks = Lists.newArrayList();

        Queue<GappedAssembly> secondaryAssemblyQueue = new ConcurrentLinkedQueue<>();
        secondaryAssemblyQueue.addAll(secondaryAssemblies);

        int secondaryAssemblyCount = secondaryAssemblyQueue.size();

        for(int i = 0; i < taskCount; ++i)
        {
            Aligner aligner = new Aligner(config, secondaryAssemblyQueue, new File(config.RefGenomeImageFile));
            alignerTasks.add(aligner);
            threadTasks.add(aligner);
        }

        SV_LOGGER.debug("splitting {} secondary assemblies across {} threads", secondaryAssemblyCount, taskCount);

        return alignerTasks;
    }

    public Aligner(final SvConfig config, final Queue<GappedAssembly> secondaryAssemblyQueue, final File index)
    {
        super("Aligner");

        mConfig = config;
        mSecondaryAssemblyQueue = secondaryAssemblyQueue;
        mSecondaryAssemblyCount = secondaryAssemblyQueue.size();

        mIndex = new BwaMemIndex(index.getAbsolutePath());
        mAligners = ThreadLocal.withInitial(() ->
        {
            final BwaMemAligner aligner = new BwaMemAligner(mIndex);
            aligner.setClip3PenaltyOption(0);
            aligner.setClip5PenaltyOption(0);
            aligner.setChunkSizeOption(10000000);
            return aligner;
        });

        mHomologySlider = new HomologySlider(mConfig.RefGenome);
        mAlignedAssemblies = Lists.newArrayList();
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mSecondaryAssemblyQueue.size();
                int processedCount = mSecondaryAssemblyCount - remainingCount;

                GappedAssembly secondaryAssembly = mSecondaryAssemblyQueue.remove();

                mPerfCounter.start();
                alignAssembly(secondaryAssembly);

                stopCheckLog(secondaryAssembly.toString(), mConfig.PerfLogTime);

                if(processedCount > 0 && (processedCount % TASK_LOG_COUNT) == 0)
                {
                    SV_LOGGER.info("aligned {} second assemblies, remaining({})", processedCount, remainingCount);
                }
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public static List<AlignedAssembly> mergeAlignedAssemblies(final List<Aligner> alignerTasks)
    {
        List<AlignedAssembly> combined = Lists.newArrayList();
        alignerTasks.forEach(x -> combined.addAll(x.mAlignedAssemblies));
        return combined;
    }

    public void alignAssembly(final GappedAssembly assembly)
    {
        AlignedAssembly alignedAssembly = align(assembly);

        // Left sliding (we will find mid-points after calling + de-duping)
        AlignedAssembly adjustedAssembly = mHomologySlider.slideHomology(alignedAssembly);

        mAlignedAssemblies.add(adjustedAssembly);
    }

    public AlignedAssembly align(final GappedAssembly assembly)
    {
        final List<Alignment> alignments = new ArrayList<>();
        int currentOffset = 0;
        for(ExtendedAssembly source : assembly.Sources)
        {
            final List<Alignment> sourceAlignments = align(source);
            for(Alignment alignment : sourceAlignments)
            {
                final int adjustedStartPosition = alignment.SequenceStartPosition + currentOffset;

                alignments.add(new Alignment(alignment.Chromosome, alignment.ReferenceStartPosition,
                        adjustedStartPosition, alignment.Length,
                        alignment.Inverted, alignment.Quality));
            }
            alignments.add(Alignment.gap(currentOffset + source.getLength()));

            currentOffset += source.getLength() + 1; // Add 1 to account for 'X's inserted in gapped assemblies.
        }
        alignments.remove(alignments.size() - 1);

        final AlignedAssembly aligned = AlignmentPostProcessor.flipIfRequired(new AlignedAssembly(assembly, alignments));
        return aligned;
    }

    public List<Alignment> align(final Sequence sequence)
    {
        final byte[] assemblyBases = sequence.getBases();
        final List<BwaMemAlignment> alignmentSet = mAligners.get().alignSeqs(List.of(assemblyBases)).get(0);

        @Nullable
        final BwaMemAlignment best = alignmentSet.stream()
                .filter(alignment -> alignment.getRefStart() != -1)
                .max(Comparator.comparingInt(BwaMemAlignment::getMapQual)
                        .thenComparingInt(BwaMemAlignment::getAlignerScore)
                        .thenComparing(Comparator.comparingInt(BwaMemAlignment::getNMismatches).reversed()))
                .orElse(null);
        if(best == null)
            return List.of(Alignment.unmapped(sequence.getLength()));

        if(SAMFlag.getFlags(best.getSamFlag()).contains(SAMFlag.READ_UNMAPPED))
            return List.of(Alignment.unmapped(sequence.getLength()));

        // Extend alignment
        List<Alignment> leftExtended = toAlignments(best);
        while(leftExtended.size() > 1 && leftExtended.get(0).isUnmapped())
        {
            final List<Alignment> newAlignment = extendAlignmentLeft(sequence, leftExtended);
            if(newAlignment == leftExtended)
                break;
            leftExtended = newAlignment;
        }

        List<Alignment> rightExtended = leftExtended;
        while(rightExtended.size() > 1 && rightExtended.get(rightExtended.size() - 1).isUnmapped())
        {
            final List<Alignment> newAlignment = extendAlignmentRight(sequence, rightExtended);
            if(newAlignment == rightExtended)
                break;
            rightExtended = newAlignment;
        }

        return cleanAlignment(unmapRubbish(fillGaps(sequence, rightExtended)));
    }

    private List<Alignment> toAlignments(final BwaMemAlignment alignment)
    {
        return toAlignments(alignment, 0);
    }

    private List<Alignment> toAlignments(final BwaMemAlignment alignment, final int sequenceOffset)
    {
        final Set<SAMFlag> flags = SAMFlag.getFlags(alignment.getSamFlag());
        final boolean inverted = flags.contains(SAMFlag.READ_REVERSE_STRAND);

        final Cigar cigar = CigarUtils.cigarFromStr(alignment.getCigar());

        final String chromosome = mIndex.getReferenceContigNames().get(alignment.getRefId());

        int referencePosition = !inverted ? alignment.getRefStart() + 1 : alignment.getRefEnd() + 1;
        int readPosition = 1 + sequenceOffset;
        final List<Alignment> alignments = new ArrayList<>();
        if(inverted)
        {
            final List<CigarElement> cigarElements = cigar.getCigarElements();
            for(int i = cigarElements.size() - 1; i >= 0; i--)
            {
                final CigarElement element = cigarElements.get(i);
                if(element.getOperator().consumesReferenceBases())
                    referencePosition -= element.getLength();
                switch(element.getOperator())
                {
                    case M:
                        alignments.add(new Alignment(chromosome, referencePosition,
                                readPosition, element.getLength(), inverted, alignment.getMapQual()));
                        break;
                    case I:
                        alignments.add(Alignment.insert(readPosition, element.getLength()));
                        break;
                    case S:
                        alignments.add(Alignment.unmapped(readPosition, element.getLength()));
                        break;
                }
                if(element.getOperator().consumesReadBases())
                    readPosition += element.getLength();
            }
        }
        else
        {
            for(CigarElement element : cigar.getCigarElements())
            {
                switch(element.getOperator())
                {
                    case M:
                        alignments.add(new Alignment(chromosome, referencePosition, readPosition, element.getLength(), inverted, alignment.getMapQual()));
                        break;
                    case I:
                        alignments.add(Alignment.insert(readPosition, element.getLength()));
                        break;
                    case S:
                        alignments.add(Alignment.unmapped(readPosition, element.getLength()));
                        break;
                }
                if(element.getOperator().consumesReadBases())
                    readPosition += element.getLength();
                if(element.getOperator().consumesReferenceBases())
                    referencePosition += element.getLength();
            }
        }
        return alignments;
    }

    private int nearbyDistance(final Alignment neighbour, final BwaMemAlignment maybeNearby)
    {
        final String maybeNearbyChromosome = mIndex.getReferenceContigNames().get(maybeNearby.getRefId());
        if(!neighbour.Chromosome.equals(maybeNearbyChromosome))
            return Integer.MAX_VALUE;

        final boolean isInverted = SAMFlag.getFlags(maybeNearby.getSamFlag()).contains(SAMFlag.READ_REVERSE_STRAND);
        if(neighbour.Inverted != isInverted)
            return Integer.MAX_VALUE;

        final int neighbourMidPoint = neighbour.ReferenceStartPosition + neighbour.Length / 2;
        final int nearbyMidPoint = ((maybeNearby.getRefStart() + 1) + (maybeNearby.getRefEnd() + 1)) / 2;
        final int distance = Math.abs(neighbourMidPoint - nearbyMidPoint);

        return distance > SvConstants.ALIGNER_NEARBY_MAX_DISTANCE ? Integer.MAX_VALUE : distance;
    }

    private List<Alignment> bestAlignmentExtension(final Sequence unmappedBases,
            final Alignment mappedNeighbour, final Direction direction)
    {
        final List<BwaMemAlignment> mappedAlignments = mAligners.get().alignSeqs(List.of(unmappedBases.getBases())).get(0).stream()
                .filter(alignment -> alignment.getSeqStart() != -1)
                .filter(alignment -> alignment.getAlignerScore() > SvConstants.ALIGNER_MIN_SCORE)
                .collect(Collectors.toList());

        if(mappedAlignments.isEmpty())
            return List.of();

        final int bestSequenceBoundary = direction == Direction.REVERSE
                ? mappedAlignments.stream().mapToInt(BwaMemAlignment::getSeqEnd).max().orElseThrow()
                : mappedAlignments.stream().mapToInt(BwaMemAlignment::getSeqStart).min().orElseThrow();

        final List<BwaMemAlignment> candidateAlignments = mappedAlignments.stream()
                .filter(alignment -> direction == Direction.REVERSE
                        ? alignment.getSeqEnd() >= bestSequenceBoundary - SvConstants.ALIGNER_EXTENSION_INSERT_TOLERANCE
                        : alignment.getSeqStart() <= bestSequenceBoundary + SvConstants.ALIGNER_EXTENSION_INSERT_TOLERANCE)
                .collect(Collectors.toList());

        // We select the best candidate as follows:
        // -> If an alignment could show a deletion (< 10k bases), choose the one that shows the smallest deletion
        // -> Otherwise, choose the alignment with the highest aligner score

        final BwaMemAlignment bestAlignment = candidateAlignments.stream()
                .min(Comparator.<BwaMemAlignment>comparingInt(alignment -> nearbyDistance(mappedNeighbour, alignment))
                        .thenComparingInt(alignment -> -alignment.getAlignerScore()))
                .orElseThrow();

        return toAlignments(bestAlignment);
    }

    private List<Alignment> extendAlignmentLeft(final Sequence sequence, final List<Alignment> existing)
    {
        if(existing.size() < 2 || !existing.get(0).isUnmapped() || existing.get(1).isUnmapped())
            return existing;

        final Alignment unmappedLeft = existing.get(0);
        final Alignment mappedLeft = existing.get(1);
        if(unmappedLeft.Length < SvConstants.ALIGNER_MIN_BASES)
            return existing;

        final Sequence unmapped = sequence.subsequence(unmappedLeft.SequenceStartPosition - 1, unmappedLeft.Length);
        if(isMonoBase(unmapped))
            return existing;

        final List<Alignment> newAlignmentBlocks = bestAlignmentExtension(unmapped, mappedLeft, Direction.REVERSE);
        if(newAlignmentBlocks.isEmpty())
            return existing;
        existing.stream().skip(1).forEach(newAlignmentBlocks::add);

        return newAlignmentBlocks;
    }

    private List<Alignment> extendAlignmentRight(final Sequence sequence, final List<Alignment> existing)
    {
        if(existing.size() < 2 || !existing.get(existing.size() - 1).isUnmapped()
                || existing.get(existing.size() - 2).isUnmapped())
            return existing;

        // Grab the first shit and go do stuff
        final Alignment unmappedRight = existing.get(existing.size() - 1);
        final Alignment mappedRight = existing.get(existing.size() - 2);

        final Sequence unmapped = sequence.subsequence(unmappedRight.SequenceStartPosition - 1, unmappedRight.Length);
        if(isMonoBase(unmapped))
            return existing;

        final List<Alignment> extensionBlocks = bestAlignmentExtension(unmapped, mappedRight, Direction.FORWARDS);
        if(extensionBlocks.isEmpty())
            return existing;

        final List<Alignment> newAlignmentBlocks = new ArrayList<>();
        existing.stream().limit(existing.size() - 1).forEach(newAlignmentBlocks::add);
        for(Alignment extension : extensionBlocks)
            newAlignmentBlocks.add(new Alignment(extension.Chromosome, extension.ReferenceStartPosition,
                    extension.SequenceStartPosition - 1 + unmappedRight.SequenceStartPosition, extension.Length,
                    extension.Inverted, extension.Quality));

        return newAlignmentBlocks;
    }

    private List<Alignment> fillGaps(final Sequence sequence, final List<Alignment> existing)
    {
        final List<Alignment> alignments = new ArrayList<>();
        for(Alignment alignment : existing)
            alignments.addAll(tryMap(sequence, alignment));
        return alignments;
    }

    private List<Alignment> tryMap(final Sequence sequence, final Alignment alignment)
    {
        if(alignment.isMapped() || alignment.Length < 30)
            return List.of(alignment);

        final int startIndex = alignment.SequenceStartPosition - 1;
        final int endIndexExcl = startIndex + alignment.Length;
        final byte[] blockBases = Arrays.copyOfRange(sequence.getBases(), startIndex, endIndexExcl);
        final List<BwaMemAlignment> candidates = mAligners.get().alignSeqs(List.of(blockBases)).get(0);
        @Nullable
        final BwaMemAlignment bestCandidate = candidates.stream()
                .max(Comparator.comparingInt(BwaMemAlignment::getAlignerScore))
                .orElse(null);
        if(bestCandidate == null || bestCandidate.getRefId() < 0)
            return List.of(alignment);

        final List<Alignment> replacements = toAlignments(bestCandidate, alignment.SequenceStartPosition - 1);
        final List<Alignment> recursiveReplacements = new ArrayList<>();

        SV_LOGGER.trace("remapping {} as {}", alignment, replacements);

        for(Alignment replacement : replacements)
        {
            if(replacement.isUnmapped() && replacement.Length >= 30)
                recursiveReplacements.addAll(tryMap(sequence, replacement));
            else
                recursiveReplacements.add(replacement);
        }
        return recursiveReplacements;
    }

    private List<Alignment> unmapRubbish(final List<Alignment> existing)
    {
        final List<Alignment> alignments = new ArrayList<>();
        for(Alignment alignment : existing)
        {
            if(alignment.isMapped() && alignment.Quality < 10)
            {
                SV_LOGGER.trace("unmapping {}", alignment);

                alignments.add(Alignment.unmapped(alignment.SequenceStartPosition, alignment.Length));
            }
            else
                alignments.add(alignment);
        }

        return alignments;
    }

    private List<Alignment> cleanAlignment(final List<Alignment> alignments)
    {
        final List<Alignment> cleaned = new ArrayList<>();
        for(int i = 0; i < alignments.size(); i++)
        {
            if(i == 0)
            {
                cleaned.add(alignments.get(i));
                continue;
            }

            final Alignment alignment = alignments.get(i);
            final Alignment previous = cleaned.get(cleaned.size() - 1);

            if(previous.isUnmapped() && alignment.isUnmapped())
                cleaned.set(cleaned.size() - 1, new Alignment(previous.Chromosome, 0, previous.SequenceStartPosition,
                        previous.Length + alignment.Length, false, 0));
            else
                cleaned.add(alignment);
        }
        return cleaned;
    }

    @Override
    public void close()
    {
        mIndex.close();
    }

    private static boolean isMonoBase(final Sequence sequence)
    {
        final float[] frequencies = baseFrequencies(sequence);
        for(float frequency : frequencies)
            if(frequency > 0.9)
                return true;
        return false;
    }

    private static float[] baseFrequencies(final Sequence sequence)
    {
        final int[] counts = baseCounts(sequence);
        int total = 0;
        for(int count : counts)
            total += count;

        final float[] frequencies = new float[counts.length];
        for(int i = 0; i < counts.length; i++)
            frequencies[i] = (float) counts[i] / total;

        return frequencies;
    }

    private static int[] baseCounts(final Sequence sequence)
    {
        final int[] counts = new int[5];
        for(byte base : sequence.getBases())
        {
            if(base == 'A')
                counts[0]++;
            else if(base == 'T')
                counts[1]++;
            else if(base == 'C')
                counts[2]++;
            else if(base == 'G')
                counts[3]++;
            else
                counts[4]++;
        }

        return counts;
    }
}
