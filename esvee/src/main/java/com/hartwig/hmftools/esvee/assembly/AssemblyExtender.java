package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.ASSEMBLY_EXTENSION_MIN_SUPPORT;
import static com.hartwig.hmftools.esvee.SvConstants.ASSEMBLY_EXTENSION_MIN_SUPPORT_FINAL;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.html.DiagramSet;
import com.hartwig.hmftools.esvee.processor.ThreadTask;
import com.hartwig.hmftools.esvee.sequence.ExtendedAssembly;
import com.hartwig.hmftools.esvee.sequence.PrimaryAssembly;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.sequence.ReadSupport;
import com.hartwig.hmftools.esvee.sequence.SupportedAssembly;
import com.hartwig.hmftools.esvee.util.Counter;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class AssemblyExtender extends ThreadTask
{
    private final Map<String,List<JunctionGroup>> mJunctionGroupMap;
    private final Queue<PrimaryAssembly> mPrimaryAssemblyQueue;
    private final int mPrimaryAssemblyCount;
    private final SvConfig mConfig;

    private final List<ExtendedAssembly> mExtendedAssemblies;

    private final SupportChecker mSupportChecker;
    private final NodeFolder mNodeFolder;
    private final AssemblyExtenderCounters mCounters;
    private final boolean mCreateDiagrams;

    private int mNextAssemblyNumber;

    public static List<AssemblyExtender> createThreadTasks(
            final Map<String,List<JunctionGroup>> junctionGroupMap, final List<PrimaryAssembly> primaryAssemblies,
            final SvConfig config, final int taskCount, final List<Thread> threadTasks)
    {
        List<AssemblyExtender> assemblyExtenderTasks = com.google.common.collect.Lists.newArrayList();

        Queue<PrimaryAssembly> primaryAssemblyQueue = new ConcurrentLinkedQueue<>();
        primaryAssemblyQueue.addAll(primaryAssemblies);

        int primaryAssemblyCount = primaryAssemblies.size();

        for(int i = 0; i < taskCount; ++i)
        {
            AssemblyExtender assemblyExtender = new AssemblyExtender(config, junctionGroupMap, primaryAssemblyQueue);
            assemblyExtenderTasks.add(assemblyExtender);
            threadTasks.add(assemblyExtender);
        }

        SV_LOGGER.debug("splitting {} primary assemblies across {} threads", primaryAssemblyCount, taskCount);

        return assemblyExtenderTasks;
    }

    public AssemblyExtender(
            final SvConfig config, final Map<String,List<JunctionGroup>> junctionGroupMap, final Queue<PrimaryAssembly> primaryAssemblyQueue)
    {
        super("AssemblyExtender");
        mJunctionGroupMap = junctionGroupMap;
        mPrimaryAssemblyQueue = primaryAssemblyQueue;
        mPrimaryAssemblyCount = primaryAssemblyQueue.size();
        mConfig = config;
        mCreateDiagrams = config.writeHtmlFiles() && config.PlotDiagrams;
        mSupportChecker = new SupportChecker();
        mCounters = new AssemblyExtenderCounters();
        mNodeFolder = new NodeFolder();
        mExtendedAssemblies = Lists.newArrayList();
        mNextAssemblyNumber = 1;

        start();
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mPrimaryAssemblyQueue.size();
                int processedCount = mPrimaryAssemblyCount - remainingCount;

                PrimaryAssembly primaryAssembly = mPrimaryAssemblyQueue.remove();

                mPerfCounter.start();
                processPrimaryAssembly(primaryAssembly);
                mPerfCounter.stop();

                if(processedCount > 0 && (processedCount % 100) == 0)
                {
                    SV_LOGGER.info("processed {} junction groups, remaining({})", processedCount, remainingCount);
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

    public List<ExtendedAssembly> extendedAssemblies() { return mExtendedAssemblies; }

    public void processPrimaryAssembly(final PrimaryAssembly assembly)
    {
        JunctionGroup junctionGroup = findJunctionGroup(assembly.OriginalJunction);

        List<Read> discordantReads = new DiscordantPairFinder().findDiscordantReads(assembly.supportingReads(), junctionGroup.candidateReads());

        mCounters.DiscordantReadsFound.add(discordantReads.size());

        List<ExtendedAssembly> leftExtended = extendLeft(assembly, discordantReads);

        leftExtended = limitBasedOnSupport(AssemblyFiltering.trimAndDeduplicate(
                mSupportChecker, leftExtended), ASSEMBLY_EXTENSION_MIN_SUPPORT);

        List<ExtendedAssembly> rightExtended = leftExtended.stream()
                .flatMap(x -> extendRight(x, discordantReads).stream()).collect(Collectors.toList());

        List<ExtendedAssembly> extendedAssemblies = limitBasedOnSupport(AssemblyFiltering.trimAndDeduplicate(
                mSupportChecker, rightExtended), ASSEMBLY_EXTENSION_MIN_SUPPORT_FINAL);

        mExtendedAssemblies.addAll(extendedAssemblies);
    }

    private JunctionGroup findJunctionGroup(final Junction junction)
    {
        List<JunctionGroup> junctionGroups = mJunctionGroupMap.get(junction.Chromosome);

        if(junctionGroups == null)
            return null;

        return junctionGroups.stream().filter(x -> x.junctions().contains(junction)).findFirst().orElse(null);
    }

    private List<ExtendedAssembly> limitBasedOnSupport(final List<ExtendedAssembly> assemblies, final int limit)
    {
        // If we've ended up with too many assemblies from this junction, raise the minimum support requirements to resolve.
        while(assemblies.size() > limit)
        {
            final int minSupport = assemblies.stream()
                    .mapToInt(candidate -> candidate.getSupportReadNames().size())
                    .min()
                    .orElseThrow();

            assemblies.removeIf(candidate -> candidate.getSupportReadNames().size() <= minSupport);
        }
        return assemblies;
    }

    private String nextAssemblyName(final String assemblyName)
    {
        return String.format("%sE%s", assemblyName, mNextAssemblyNumber++);
    }

    /*
    private List<Pair<Read, Read>> findMates(final Collection<Read> reads)
    {
        return mSAMSource.streamMates(reads).collect(Collectors.toList());

        as per the SAMSource interface:
        final Stream<RegionOfInterest> mateRegionsUnmerged = records.stream()
                .map(record -> record.isMateMapped()
                        ? new RegionOfInterest(record.getMateChromosome(), record.getMateAlignmentStart(),
                        record.getMateAlignmentStart() + record.getLength())
                        : new RegionOfInterest(record.getChromosome(), record.getAlignmentStart(),
                                record.getAlignmentStart() + record.getLength()));
        final List<RegionOfInterest> mappedMateRegions = RegionOfInterest.tryMerge(mateRegionsUnmerged::iterator);

        final Map<String, List<Record>> recordsByName = records.stream()
                .collect(MultiMapCollector.keyed(Record::getName));

        return mappedMateRegions.stream()
                .map(mateRegion -> streamReadsContaining(mateRegion.Chromosome, mateRegion.Start, mateRegion.End)
                        .flatMap(record -> recordsByName.getOrDefault(record.getName(), List.of()).stream()
                                .filter(match -> match.isFirstOfPair() != record.isFirstOfPair())
                                .map(match -> Pair.of(match, record))))
                .reduce(Stream::concat).orElse(Stream.of());
    }
    */

    private List<Read> findMateReads(final SupportedAssembly assembly, byte requiredOrientation)
    {
        List<Read> mateReads = Lists.newArrayList();

        for(Read read : assembly.supportingReads())
        {
            if(!read.hasMateSet())
                continue;

            if(read.orientation() != requiredOrientation)
                continue;

            Read mateRead = read.mateRead();

            if(assembly.containsSupport(mateRead)) // a little inefficient?
                continue;

            // CHECK: flip read if pair is an INV

            mateReads.add(mateRead);
        }

        if(requiredOrientation == NEG_ORIENT)
            Collections.sort(mateReads, Comparator.comparingInt(Read::getUnclippedEnd).reversed());
        else
            Collections.sort(mateReads, Comparator.comparingInt(Read::getUnclippedStart));

        // CHECK: was .sorted(Comparator.comparingInt(pair -> assembly.getSupportIndex(pair.getLeft())))

        return mateReads;
    }

    private List<ExtendedAssembly> extendLeft(final SupportedAssembly assembly, final List<Read> discordantReads)
    {
        List<Read> mateReads = findMateReads(assembly, NEG_ORIENT);

        /*
        List<Pair<Read, Read>> pairedMates = findMates(assembly.getSupportRecords().stream()
                .filter(x -> x.negativeStrand())
                .collect(Collectors.toList()));

        // note that this method flips the right read if this pair is an INV
        final List<Read> mates = pairedMates.stream()
                .filter(pair -> !assembly.containsSupport(pair.getRight()))
                .filter(pair -> AlignmentFilters.isRecordAverageQualityAbove(pair.getRight().getBaseQuality(), SvConstants.AVG_BASE_QUAL_THRESHOLD))
                .map(pair -> pair.getLeft().positiveStrand() == pair.getRight().positiveStrand()
                        ? ReadUtils.flipRead(pair.getRight())
                        : pair.getRight())
                .sorted(Comparator.comparingInt(Read::getUnclippedEnd).reversed())
                .collect(Collectors.toList());
        */

        mCounters.LeftMates.add(mateReads.size());

        discordantReads.sort(Comparator.comparingInt(Read::getUnclippedEnd).reversed());

        List<Read> applicableDiscordantReads = discordantReads.stream().filter(x -> !x.negativeStrand()).collect(Collectors.toList());

        final Map<Read, Integer> checkStartIndices = new LinkedHashMap<>();
        for(Read read : Iterables.concat(mateReads, applicableDiscordantReads))
        {
            // Get the support index of this record or its mate
            List<ReadSupport> readSupports = assembly.getReadSupport(read.getName());

            Integer existingSupportIndex = null;

            // CHECK: take the lowest index including read length - why?
            if(readSupports != null)
            {
                for(ReadSupport readSupport : readSupports)
                {
                    int combinedLength = readSupport.Index + readSupport.Read.getLength();

                    if(existingSupportIndex == null || combinedLength < existingSupportIndex)
                        existingSupportIndex = combinedLength;
                }

            }

            /*
            Integer existingSupportIndex = assembly.getSupport(read.getName()).stream()
                    .map(entry -> entry.getValue() + entry.getKey().getLength())
                    .min(Comparator.naturalOrder())
                    .orElse(null);
            */

            int minDepth = existingSupportIndex == null
                    ? assembly.getLength() - read.getLength()
                    : Math.max(assembly.getLength() - read.getLength(), assembly.getLength() - existingSupportIndex);
            checkStartIndices.put(read, minDepth);
        }

        Direction assemblyReadDirection = Direction.REVERSE;

        Direction mateReadDirection = Direction.REVERSE;

        List<ExtendedAssembly> extended = extendAssembly(assembly, assemblyReadDirection, mateReads, discordantReads,
                checkStartIndices, mateReadDirection, mCounters.LeftMatesAssembled, mCounters.LeftDiscordantReadsAssembled);

        mCounters.ExtendLeftAssemblies.add(extended.size());

        return extended;
    }

    private List<ExtendedAssembly> extendRight(final SupportedAssembly assembly, final List<Read> discordantReads)
    {
        /*
        final List<Pair<Read, Read>> pairedMates = findMates(assembly.getSupportRecords().stream()
                .filter(record -> !record.negativeStrand())
                .collect(Collectors.toList()));

        // CHECK: reason for flipping?
        // CHECK: why checking avg base qual again when already checked for all reads?
        final List<Read> mates = pairedMates.stream()
                .filter(pair -> !assembly.containsSupport(pair.getRight()))
                .filter(pair -> AlignmentFilters.isRecordAverageQualityAbove(pair.getRight().getBaseQuality(), SvConstants.AVG_BASE_QUAL_THRESHOLD))
                .sorted(Comparator.comparingInt(pair -> assembly.getSupportIndex(pair.getLeft())))
                .map(pair -> pair.getLeft().positiveStrand() == pair.getRight().positiveStrand()
                        ? ReadUtils.flipRead(pair.getRight()) : pair.getRight())
                .collect(Collectors.toList());
        */

        List<Read> mateReads = findMateReads(assembly, POS_ORIENT);

        mCounters.RightMates.add(mateReads.size());

        discordantReads.sort(Comparator.comparingInt(Read::getUnclippedStart));

        List<Read> applicableDiscordantReads = discordantReads.stream().filter(x -> x.negativeStrand()).collect(Collectors.toList());

        final Map<Read, Integer> checkStartIndices = new LinkedHashMap<>();

        for(Read read : Iterables.concat(mateReads, applicableDiscordantReads))
        {
            // Get the support index of this record or its mate
            Integer existingSupportIndex = null;

            List<ReadSupport> readSupports = assembly.getReadSupport(read.getName());

            if(readSupports != null)
            {
                // CHECK: shouldn't be null?
                existingSupportIndex = readSupports.stream().mapToInt(x -> x.Index).max().orElse(-1);
            }

            /*
            Integer existingSupportIndex = assembly.getSupport(read.getName()).stream()
                    .map(Map.Entry::getValue)
                    .max(Comparator.naturalOrder())
                    .orElse(null);
             */

            int minDepth = existingSupportIndex == null
                    ? assembly.getLength() - read.getLength()
                    : Math.max(assembly.getLength() - read.getLength(), existingSupportIndex);

            checkStartIndices.put(read, minDepth);
        }

        Direction assemblyReadDirection = Direction.FORWARDS;

        Direction mateReadDirection = Direction.FORWARDS;

        List<ExtendedAssembly> extended = extendAssembly(
                assembly, assemblyReadDirection, mateReads, applicableDiscordantReads,
                checkStartIndices, mateReadDirection, mCounters.RightMatesAssembled, mCounters.RightDiscordantReadsAssembled);

        mCounters.ExtendRightAssemblies.add(extended.size());

        return extended;
    }

    private List<ExtendedAssembly> extendAssembly(
            final SupportedAssembly assembly, final Direction assemblyDirection,
            final List<Read> mateAlignments,
            final List<Read> discordantAlignments,
            final Map<Read, Integer> alignmentMinDepth,
            final Direction alignmentDirection,
            final Counter mateAttachCounter,
            final Counter discordantAttachCounter)
    {
        if(mateAlignments.isEmpty() && discordantAlignments.isEmpty())
        {
            if(assembly instanceof ExtendedAssembly)
                return List.of((ExtendedAssembly) assembly);

            ExtendedAssembly newAssembly = new ExtendedAssembly(assembly.Name, assembly.Assembly, assembly);

            assembly.readSupport().forEach(x -> newAssembly.addEvidenceAt(x.Read, x.Index));
            return List.of(newAssembly);
        }

        final HeadNode existing = HeadNode.create(assembly, assemblyDirection);

        final List<Read> potentialNewSupport = new ArrayList<>();
        extendAssembly(existing, potentialNewSupport, mateAlignments, alignmentMinDepth, alignmentDirection, mateAttachCounter);
        extendAssembly(existing, potentialNewSupport, discordantAlignments, alignmentMinDepth, alignmentDirection, discordantAttachCounter);

        final Map<Read, Set<Integer>> supportStartIndices = potentialSupportIndices(existing);

        final String diagramSetName = "Extension " + assemblyDirection.name().toLowerCase();
        @Nullable
        final DiagramSet diagrams = simplifyGraph(diagramSetName, existing);

        return existing.flatten().stream()
                .map(newAssembly ->
                {
                    final String directionCorrected = assemblyDirection == Direction.FORWARDS
                            ? newAssembly
                            : new StringBuilder(newAssembly).reverse().toString();

                    final ExtendedAssembly newCandidate = new ExtendedAssembly(nextAssemblyName(assembly.Name), directionCorrected, assembly);
                    newCandidate.addDiagrams(diagrams);

                    reAddSupport(newCandidate, assembly, potentialNewSupport, supportStartIndices,
                            assemblyDirection == Direction.FORWARDS);

                    return newCandidate;
                }).collect(Collectors.toList());
    }

    private void extendAssembly(
            final HeadNode existing, final List<Read> potentialNewSupport, final List<Read> alignments,
            final Map<Read, Integer> alignmentMinDepth, final Direction alignmentDirection, final Counter attachCounter)
    {
        for(Read read : alignments)
        {
            @Nullable
            final Integer depth = alignmentMinDepth.get(read);
            final HeadNode toAttach = alignmentAsGraph(read, alignmentDirection);

            if(existing.attach(toAttach, SvConstants.ASSEMBLY_EXTENSION_MIN_MATCH_BASES, SvConstants.ASSEMBLY_EXTENSION_MAX_MISMATCH,
                    Objects.requireNonNullElse(depth, 0), Integer.MAX_VALUE))
            {
                potentialNewSupport.add(read);
                attachCounter.add(1);
            }
        }
    }

    // Once we've built the graph, try and work out where in the flattened output different records may appear
    private Map<Read, Set<Integer>> potentialSupportIndices(final HeadNode head)
    {
        final Map<Read, Set<Integer>> candidateLocations = new HashMap<>();
        final Deque<Pair<Node, Integer>> worklist = new ArrayDeque<>();
        final Set<Pair<Node, Integer>> seen = new HashSet<>();
        for(Node successor : head.successors())
            worklist.add(Pair.of(successor, 0));

        while(!worklist.isEmpty())
        {
            final Pair<Node, Integer> pair = worklist.poll();
            if(!seen.add(pair))
                continue;

            final Node node = pair.getLeft();
            final int depth = pair.getRight();

            for(Node.Support support : node.Support)
                if(support.ReadIndex <= 2)
                    candidateLocations.computeIfAbsent(support.Read, r -> new LinkedHashSet<>()).add(depth - support.ReadIndex);

            for(Node successor : node.successors())
                worklist.add(Pair.of(successor, depth + 1));
        }

        return candidateLocations;
    }

    private void reAddSupport(
            final ExtendedAssembly assembly, final SupportedAssembly original, final List<Read> potentialSupport,
            final Map<Read, Set<Integer>> supportStartIndices, final boolean isForwards)
    {
        if(potentialSupport.size() > 1_000)
            return;

        final int firstIndex = assembly.Assembly.indexOf(original.Assembly);
        final int lastIndex = assembly.Assembly.indexOf(original.Assembly, firstIndex);
        if(firstIndex == -1 || firstIndex != lastIndex)
            addSupportNoOffset(assembly, original);
        else
            addSupportWithOffset(assembly, original, firstIndex);

        for(Read read : potentialSupport)
        {
            if(assembly.containsSupport(read))
                continue;

            @Nullable
            final Set<Integer> candidates = supportStartIndices.get(read);

            if(candidates == null)
            {
                assembly.tryAddSupport(mSupportChecker, read);
            }
            else
            {
                for(int rawIndex : candidates)
                {
                    final int checkIndex = isForwards
                            ? rawIndex
                            : assembly.getLength() - rawIndex - read.getLength();
                    if(mSupportChecker.WeakSupport.supportsAt(assembly, read, checkIndex))
                    {
                        assembly.addEvidenceAt(read, checkIndex);
                        break;
                    }
                }
            }
        }
    }

    private void addSupportWithOffset(final ExtendedAssembly assembly, final SupportedAssembly original, final int offset)
    {
        original.readSupport().forEach(x ->
        {
            if(mSupportChecker.WeakSupport.supportsAt(assembly, x.Read, x.Index + offset))
            {
                assembly.addEvidenceAt(x.Read, x.Index + offset);
            }
            else
            {
                mSupportChecker.supportsAtIndex(original, x.Read, 20, 2, x.Index);
                mSupportChecker.supportsAtIndex(original, x.Read, 20, 2, -32);
            }
        });
    }

    private void addSupportNoOffset(final ExtendedAssembly assembly, final SupportedAssembly original)
    {
        original.readSupport().forEach(x -> assembly.tryAddSupport(mSupportChecker, x.Read));
    }

    @Nullable
    private DiagramSet simplifyGraph(final String diagramSetName, final HeadNode node)
    {
        final DiagramSet diagrams;
        if(mCreateDiagrams)
        {
            diagrams = new DiagramSet(diagramSetName);
            diagrams.add("Attachment", node.toDiagram());
            mNodeFolder.prepruneNodes(node);
            diagrams.add("Pre-folding Prune", node.toDiagram());
            mNodeFolder.foldPaths(node);
            diagrams.add("Folding", node.toDiagram());
            node.pruneNodes();
            diagrams.add("Pruning", node.toDiagram());
        }
        else
        {
            diagrams = null;
            mNodeFolder.foldPaths(node);
            node.pruneNodes();
        }

        return diagrams;
    }

    private HeadNode alignmentAsGraph(final Read alignment, final Direction orientation)
    {
        final HeadNode root = new HeadNode();
        Node current = root;
        for(int i = 0; i < alignment.getLength(); i++)
        {
            final int index = orientation == Direction.FORWARDS ? i : alignment.getLength() - 1 - i;
            final byte base = alignment.getBases()[index];
            final int quality = alignment.getBaseQuality()[index];

            final Node newNode = new Node((char) base);
            newNode.MaxQuality = newNode.Quality = quality;
            newNode.Support = new ArrayList<>();
            newNode.Support.add(new Node.Support(alignment, i));
            current.setNext(newNode);
            current = newNode;
        }

        return root;
    }
}
