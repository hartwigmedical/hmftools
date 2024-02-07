package com.hartwig.hmftools.esvee.old;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_WEAK_SUPPORT_MIN_BASES;
import static com.hartwig.hmftools.esvee.SvConstants.READ_FILTER_MIN_ALIGNED_BASES;
import static com.hartwig.hmftools.esvee.old.ReadFilters.hasAcceptableMapQ;
import static com.hartwig.hmftools.esvee.old.ReadFilters.isRecordAverageQualityPastJunctionAbove;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.read.ReadFilters;
import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.read.Read;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class PrimaryAssembler
{
    private final SvConfig mConfig;

    private final Junction mJunction;
    
    private final SupportChecker mSupportChecker;
    private final NodeFolder mNodeFolder;
    private final PrimaryAssemblerCounters mCounters;

    private int mNextAssemblyNumber = 1;

    public PrimaryAssembler(final SvConfig config, final Junction junction)
    {
        mConfig = config;
        mSupportChecker = new SupportChecker();
        mNodeFolder = new NodeFolder();
        mCounters = new PrimaryAssemblerCounters();
        mJunction = junction;
    }

    public PrimaryAssemblerCounters getCounters()
    {
        return mCounters;
    }

    public List<PrimaryAssembly> processJunction(final List<Read> rawReads)
    {
        final List<Read> realignedReads = rawReads.stream()
                .map(alignment -> realignForJunction(alignment, mJunction))
                .collect(Collectors.toList());

        final List<Read> withLowQAlignments = realignedReads.stream()
                .filter(alignment -> ReadFilters.recordSoftClipsNearJunction(alignment, mJunction)) // mCounters.ReadsSoftClippedAtJunction
                .collect(Collectors.toList());

        final List<Read> filteredAlignments = withLowQAlignments.stream()
                .filter(alignment -> isRecordAverageQualityPastJunctionAbove(alignment, mJunction, SvConstants.AVG_BASE_QUAL_THRESHOLD)) // mCounters.ReadsPassingJunctionQualityThreshold
                .filter(alignment -> hasAcceptableMapQ(alignment, SvConstants.READ_FILTER_MIN_JUNCTION_MAPQ)) // mCounters.HasAcceptableMapQ
                .filter(x -> isNotBadlyMapped(x)) // mCounters.WellMapped
                .collect(Collectors.toList());

        if(filteredAlignments.isEmpty())
            return List.of(); // There are no reads of acceptable quality supporting this junction

        final List<PrimaryAssembly> initialAssemblies = createInitialAssemblies(filteredAlignments);

        // mCounters.InitialAssemblies.add(initialAssemblies.size());

        final List<PrimaryAssembly> extendedInitial = SvConstants.EXTEND_PRIMARIES
                 ? extendInitial(withLowQAlignments, initialAssemblies) : initialAssemblies;

        final List<PrimaryAssembly> dedupedInitial = AssemblyFiltering.trimAndDeduplicate(mSupportChecker, extendedInitial);
        // mCounters.DedupedInitialAssemblies.add(dedupedInitial.size());
        
        final List<PrimaryAssembly> anchored = createAnchors(realignedReads, dedupedInitial);

        for(PrimaryAssembly assembly : anchored)
        {
            for(Read alignment : withLowQAlignments)
            {
                if(!assembly.containsSupport(alignment))
                {
                    @Nullable
                    final Integer supportIndex = mSupportChecker.WeakSupport.bestSupportIndex(assembly, (Sequence)alignment, 50);
                    if(supportIndex != null)
                        assembly.addEvidenceAt(alignment, supportIndex);
                }
            }
        }

        final List<PrimaryAssembly> assemblies = AssemblyFiltering.trimAndDeduplicate(mSupportChecker, anchored);
        // mCounters.DedupedAnchoredAssemblies.add(assemblies.size());

        // final JunctionMetrics junctionMetrics = new JunctionMetrics(mJunction.Chromosome, mJunction.Position, mJunction.direction(), mCounters);
        // assemblies.forEach(assembly -> assembly.addErrata(junctionMetrics));

        return assemblies;
    }

    private static boolean isNotBadlyMapped(final Read read)
    {
        return !isBadlyMapped(read);
    }

    private static boolean isBadlyMapped(final Read read)
    {
        if(!isDiscordant(read))
            return false;

        final int mappedSize = read.cigarElements().stream()
                .filter(element -> element.getOperator() == CigarOperator.M)
                .mapToInt(CigarElement::getLength)
                .sum();

        int indelCount = read.cigarElements().stream()
                .filter(element -> element.getOperator() == CigarOperator.D || element.getOperator() == CigarOperator.I)
                .mapToInt(CigarElement::getLength)
                .sum();

        // CHECK
        int nmCount = max(read.numberOfEvents(), indelCount);

        final int mismatchedBases = nmCount - indelCount;
        final int matchedBases = mappedSize - mismatchedBases;

        if(matchedBases < READ_FILTER_MIN_ALIGNED_BASES)
            return true;

        if(indelCount > 0 && mismatchedBases >= 3)
            return true;


        try
        {
            final float[] stats = mappedBaseStats(read);

            int countAbove70 = 0;
            int countAbove35 = 0;
            for(float frequency : stats)
            {
                if(frequency > 0.7f)
                    countAbove70++;
                if(frequency > 0.35f)
                    countAbove35++;
            }

            if(countAbove70 >= 1 || countAbove35 >= 2)
                return true;

            return false;
        }
        catch(Exception e)
        {
            SV_LOGGER.error("failed to handle read mapping stats: {}", read.toString());
            return true;
        }
    }

    private static float[] mappedBaseStats(final Read read)
    {
        final int[] baseCount = new int[5];

        int readPosition = 1;
        for(CigarElement element : read.cigarElements())
        {
            if(element.getOperator() != CigarOperator.M)
            {
                if(element.getOperator().consumesReadBases())
                    readPosition += element.getLength();
                continue;
            }

            for(int i = 0; i < element.getLength(); i++)
            {
                byte base = read.getBases()[readPosition + i - 1];
                baseCount[baseToIndex(base)]++;
            }
        }

        final float totalBases = Arrays.stream(baseCount).sum();
        final float[] baseFrequency = new float[baseCount.length];
        for(int i = 0; i < baseCount.length; i++)
        {
            baseFrequency[i] = baseCount[i] / totalBases;
        }

        return baseFrequency;
    }

    private static int baseToIndex(final byte base)
    {
        switch(base)
        {
            case 'A':
                return 0;
            case 'T':
                return 1;
            case 'C':
                return 2;
            case 'G':
                return 3;
            default:
                return 4;
        }
    }

    private List<PrimaryAssembly> extendInitial(final List<Read> alignments, final List<PrimaryAssembly> assemblies)
    {
        return assemblies.stream()
                .map(assembly -> extendInitial(alignments, assembly, mJunction.direction()))
                .collect(Collectors.toList());
    }

    /** There may be alignments that can extend the assembly but that are too noisy to be used during initial construction.
     * Examples of these types of alignments may be, for example, ones with larger soft-clips that have resulted in unacceptably low MapQ.
     * Extension in this manner is not supposed to create new candidates, so we will always choose the "best" result after pruning. */
    private PrimaryAssembly extendInitial(final List<Read> reads, final PrimaryAssembly assembly, final Direction direction)
    {
        HeadNode graph = HeadNode.create(assembly, direction);

        List<Read> supportReads = assembly.supportingReads();

        for(Read read : reads)
        {
            if(supportReads.contains(read))
                continue;

            if(!mSupportChecker.WeakSupport.supports(assembly, (Sequence)read))
                continue; // PERF: This should be supports-at

            graph = HeadNode.combine(graph, HeadNode.create(read, assembly.AnchorPosition, direction));
        }

        final List<String> flattened = graph.flatten();

        if(flattened.size() * reads.size() > 100_000)
        {
            //throw new JunctionProcessingException("Too many flattened assemblies or alignments!");
            SV_LOGGER.info("{} got {} extensions & {} alignments for a product of {}",
                    assembly.getName(), flattened.size(), reads.size(), flattened.size() * reads.size());
        }

        return Stream.concat(Stream.of(assembly), flattened.stream()
                .map(assemblyBases ->
                {
                    if(direction == Direction.REVERSE)
                        assemblyBases = new StringBuilder(assemblyBases).reverse().toString();

                    final int anchorPositionInAssembly = direction == Direction.FORWARDS ? 1 : assemblyBases.length() - 1;

                    final PrimaryAssembly newAssembly = new PrimaryAssembly(
                            nextAssemblyName(), assemblyBases, mJunction, assembly.AnchorChromosome, assembly.AnchorPosition, anchorPositionInAssembly);

                    for(Read read : reads)
                    {
                        // To support the assembly we need to either be fully contained in the assembly, or to support
                        // it with our back half if we're a forwards junction / front half if we're a backwards junction.
                        final int minSupportIndex = mJunction.direction() == Direction.FORWARDS
                                ? -read.basesLength()
                                : 0;
                        final int maxSupportIndex = mJunction.direction() == Direction.FORWARDS
                                ? min(0, newAssembly.Assembly.length() - read.basesLength())
                                : newAssembly.Assembly.length();

                        @Nullable
                        final Integer supportIndex = mSupportChecker.StrongSupport.supportIndex(newAssembly, (Sequence)read, 3, minSupportIndex, maxSupportIndex);
                        if(supportIndex != null)
                            newAssembly.addEvidenceAt(read, supportIndex);
                    }
                    return newAssembly;
                }))
                .max(Comparator.comparingInt(newAssembly -> newAssembly.getSupportReadNames().size()))
                .orElseThrow();
    }

    private String nextAssemblyName()
    {
        return String.format("%s:%s%s:%s", mJunction.Chromosome, mJunction.Position,
                mJunction.direction() == Direction.FORWARDS ? "F" : "R", mNextAssemblyNumber++);
    }

    // convert indels near the junction to soft-clips
    public static Read realignForJunction(final Read read, final Junction junction)
    {
        boolean justHadIndel = false;
        boolean wasRightNearJunction = false;
        int referencePosition = read.alignmentStart();
        for(int i = 0; i < read.getCigar().numCigarElements(); i++)
        {
            final CigarElement element = read.getCigar().getCigarElement(i);

            final boolean isIndel = element.getOperator() == CigarOperator.D || element.getOperator() == CigarOperator.I;
            final boolean isLeftNearJunction = Math.abs(referencePosition - junction.position()) <= 2;
            final boolean isRightNearJunction = Math.abs(referencePosition + element.getLength() - junction.position()) <= 2;

            if(junction.direction() == Direction.REVERSE && isLeftNearJunction && justHadIndel)
            {
                // Soft-clip everything before this cigar element (excl)
                int softClippedLength = 0;
                for(int j = 0; j < i; j++)
                    if(read.getCigar().getCigarElement(j).getOperator().consumesReadBases())
                        softClippedLength += read.getCigar().getCigarElement(j).getLength();

                final List<CigarElement> newElements = new ArrayList<>();
                newElements.add(new CigarElement(softClippedLength, CigarOperator.S));
                for(int j = i; j < read.getCigar().numCigarElements(); j++)
                    newElements.add(read.getCigar().getCigarElement(j));

                /*
                final Read newAlignment = read.copyRecord();
                newAlignment.setAlignmentStart(referencePosition);
                newAlignment.setCigar(new Cigar(newElements));
                return newAlignment;
                */

                read.bamRecord().setAlignmentStart(referencePosition);
                read.bamRecord().setCigar(new Cigar(newElements));

                return read;
            }

            if(junction.direction() == Direction.FORWARDS && wasRightNearJunction && isIndel)
            {
                // Soft-clip everything after this cigar element (incl)
                int softClippedLength = 0;
                for(int j = i; j < read.getCigar().numCigarElements(); j++)
                    if(read.getCigar().getCigarElement(j).getOperator().consumesReadBases())
                        softClippedLength += read.getCigar().getCigarElement(j).getLength();

                final List<CigarElement> newElements = new ArrayList<>();
                for(int j = 0; j < i; j++)
                    newElements.add(read.getCigar().getCigarElement(j));
                newElements.add(new CigarElement(softClippedLength, CigarOperator.S));

                read.bamRecord().setCigar(new Cigar(newElements));
                /*
                final Read newAlignment = read.copyRecord();
                newAlignment.setCigar(new Cigar(newElements));
                return newAlignment;
                */

                return read;
            }

            if(element.getOperator().consumesReferenceBases())
                referencePosition += element.getLength();
            justHadIndel = isIndel;
            wasRightNearJunction = isRightNearJunction;
        }

        return read;
    }

    private List<PrimaryAssembly> createInitialAssemblies(final List<Read> alignments)
    {
        final HeadNode combinedForwards = alignments.stream()
                .filter(alignment -> alignment.chromosome().equals(mJunction.Chromosome))
                .map(alignment -> HeadNode.create(alignment, mJunction.Position, mJunction.direction()))
                .filter(Objects::nonNull)
                .reduce(HeadNode::combine)
                .orElseThrow();

        simplifyGraph(combinedForwards, false);

        final List<String> flattened = combinedForwards.flatten();

        // mCounters.FlattenedInitial.add(flattened.size());

        final List<PrimaryAssembly> candidateAssemblies = flattened.stream()
                .filter(assemblyString -> assemblyString.length() >= 10)
                .map(assemblyString ->
                {
                    final String orientedAssembly = mJunction.direction() == Direction.FORWARDS
                            ? assemblyString
                            : new StringBuilder(assemblyString).reverse().toString();
                    final int anchorPositionInAssembly = mJunction.direction() == Direction.FORWARDS
                            ? 0
                            : orientedAssembly.length() - 1;

                    return new PrimaryAssembly(
                            nextAssemblyName(), orientedAssembly, mJunction, mJunction.Chromosome, mJunction.Position, anchorPositionInAssembly);
                })
                .collect(Collectors.toList());

        for(PrimaryAssembly assembly : candidateAssemblies)
        {
            for(Read read : alignments)
            {
                // To support the assembly we need to either be fully contained in the assembly, or to support
                // it with our back half if we're a forwards junction / front half if we're a backwards junction.
                final int minSupportIndex = mJunction.direction() == Direction.FORWARDS
                        ? -read.basesLength()
                        : 0;
                final int maxSupportIndex = mJunction.direction() == Direction.FORWARDS
                        ? min(0, assembly.Assembly.length() - read.basesLength())
                        : assembly.Assembly.length();

                @Nullable
                final Integer supportIndex = mSupportChecker.WeakSupport.supportIndex(
                        assembly, (Sequence)read, PRIMARY_ASSEMBLY_WEAK_SUPPORT_MIN_BASES, minSupportIndex, maxSupportIndex);

                if(supportIndex != null)
                    assembly.addEvidenceAt(read, supportIndex);
            }
        }

        return candidateAssemblies.stream()
                .filter(assembly -> assembly.readSupportCount() != 0)
                .collect(Collectors.toList());
    }

    private void simplifyGraph(final HeadNode node, final boolean aggressive)
    {
        mNodeFolder.foldPaths(node);

        if(aggressive)
            node.pruneNodesAggressive();
        else
            node.pruneNodes();
    }

    private List<PrimaryAssembly> createAnchors(final List<Read> alignments, final List<PrimaryAssembly> initialAssemblies)
    {
        final Map<Read, HeadNode> reverseSequences = new HashMap<>();
        for(Read read : alignments)
            reverseSequences.put(read, HeadNode.create(read, mJunction.Position, mJunction.direction().opposite()));

        final List<PrimaryAssembly> anchored = initialAssemblies.stream()
                .flatMap(candidateAssembly -> createAnchor(reverseSequences, candidateAssembly).stream())
                .collect(Collectors.toList());

        // mCounters.AnchoredAssemblies.add(anchored.size());
        return anchored;
    }

    private List<PrimaryAssembly> createAnchor(final Map<Read, HeadNode> reverseSequences, final PrimaryAssembly initialAssembly)
    {
        HeadNode anchor = null;

        for(Read support : initialAssembly.supportingReads())
        {
            HeadNode sequence = reverseSequences.get(support);

            if(anchor == null)
                anchor = sequence.deepCopy();
            else
                anchor = HeadNode.combine(anchor, sequence, false);
        }

        assert anchor != null;
        anchor.sortSupport();
        anchor = anchor.deepCopy();

        // CHECK: this used to also create Diagrams, so not completely sure it is required
        simplifyGraph(anchor, false);

        final List<String> flattened = anchor.flatten();
        // mCounters.FlattenedAnchors.add(flattened.size());
        final List<PrimaryAssembly> anchoredAssemblies = new ArrayList<>();

        for(String flattenedAssembly : flattened)
        {
            final boolean isForwards = mJunction.direction() == Direction.FORWARDS;
            final String anchoredAssembly = isForwards
                    ? new StringBuilder(flattenedAssembly.substring(1)).reverse() + initialAssembly.Assembly
                    : initialAssembly.Assembly + flattenedAssembly.substring(1);

            final int offsetSize = isForwards ? flattenedAssembly.length() - 1 : 0;
            final int anchorPositionInAssembly = initialAssembly.AnchorPositionInAssembly + offsetSize;

            // Don't use candidate to construct, as we want to re-evaluate evidence
            PrimaryAssembly assembly = new PrimaryAssembly(
                    nextAssemblyName(), anchoredAssembly,  mJunction, mJunction.Chromosome, mJunction.Position, anchorPositionInAssembly);

            for(ReadSupport readSupport : initialAssembly.readSupport())
            {
                Read read = readSupport.Read;
                int newSupportIndex;
                if(isForwards)
                    newSupportIndex = readSupport.Index + (anchoredAssembly.length() - initialAssembly.Assembly.length());
                else
                    newSupportIndex = readSupport.Index;

                if(mSupportChecker.WeakSupport.supportsAt(assembly, (Sequence)read, newSupportIndex))
                    assembly.addEvidenceAt(read, newSupportIndex);
                else
                    assembly.tryAddSupport(mSupportChecker, read);
            }

            if(assembly.getSupportReadNames().size() > SvConstants.MIN_READS_SUPPORT_ASSEMBLY)
                anchoredAssemblies.add(assembly);
        }

        return anchoredAssemblies;
    }
}
