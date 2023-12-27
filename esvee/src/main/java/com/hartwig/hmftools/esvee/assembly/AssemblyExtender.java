package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.esvee.Direction;
import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.models.DiagramSet;
import com.hartwig.hmftools.esvee.models.ExtendedAssembly;
import com.hartwig.hmftools.esvee.models.PrimaryAssembly;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.models.SupportedAssembly;
import com.hartwig.hmftools.esvee.sam.SAMSource;
import com.hartwig.hmftools.esvee.util.Counter;
import com.hartwig.hmftools.esvee.processor.Problem;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class AssemblyExtender
{
    private final SAMSource mSAMSource;
    private final SupportChecker mSupportChecker;
    private final NodeFolder mNodeFolder;
    private final PrimaryAssembly mPrimary;
    private final AssemblyExtenderCounters mCounters = new AssemblyExtenderCounters();
    private final boolean mCreateDiagrams;

    private int mNextAssemblyNumber = 1;

    @Nullable
    public static List<ExtendedAssembly> process(
            final Context context, final PrimaryAssembly assembly, final AssemblyExtenderCounters counters)
    {
        final boolean createDiagrams = context.Config.writeHtmlFiles() && context.Config.PlotDiagrams;

        final AssemblyExtender assembler = new AssemblyExtender(
                context.SAMSource, context.SupportChecker, assembly, createDiagrams);

        try
        {
            final List<ExtendedAssembly> results = assembler.extend(assembly);
            for(ExtendedAssembly extended : results)
            {
                extended.addErrata(assembly.getAllErrata());
                extended.addErrata(assembly);
                extended.addErrata(assembler.getCounters());
            }
            return results;
        }
        catch(final Throwable throwable)
        {
            final Problem problem = new Problem("Failure during extension", throwable, assembly);
            context.Problems.add(problem);

            SV_LOGGER.warn("{}", problem);

            return null;
        }
        finally
        {
            counters.add(assembler.getCounters());
        }
    }

    private AssemblyExtender(
            final SAMSource source, final SupportChecker supportChecker, final PrimaryAssembly primary, final boolean createDiagrams)
    {
        mSAMSource = source;
        mSupportChecker = supportChecker;
        mNodeFolder = new NodeFolder();
        mPrimary = primary;
        mCreateDiagrams = createDiagrams;
    }

    public AssemblyExtenderCounters getCounters()
    {
        return mCounters;
    }

    public List<ExtendedAssembly> extend(final PrimaryAssembly assembly)
    {
        return doExtend(assembly);
    }

    private List<ExtendedAssembly> limitBasedOnSupport(final List<ExtendedAssembly> assemblies, final int limit)
    {
        // If we've ended up with too many assemblies from this junction, raise the minimum support requirements to resolve.
        while(assemblies.size() > limit)
        {
            final int minSupport = assemblies.stream()
                    .mapToInt(candidate -> candidate.getSupportFragments().size())
                    .min()
                    .orElseThrow();

            assemblies.removeIf(candidate -> candidate.getSupportFragments().size() <= minSupport);
        }
        return assemblies;
    }

    private List<ExtendedAssembly> doExtend(final PrimaryAssembly assembly)
    {
        final List<Record> discordantReads = new DiscordantPairFinder(mSAMSource).findDiscordantReads(assembly.getSupportRecords());

        mCounters.DiscordantReadsFound.add(discordantReads.size());

        final List<ExtendedAssembly> leftExtended = limitBasedOnSupport(
                AssemblyFiltering.trimAndDeduplicate(mSupportChecker, extendLeft(assembly, discordantReads)), 8);

        final List<ExtendedAssembly> rightExtended = leftExtended.stream()
                .flatMap(l -> extendRight(l, discordantReads).stream())
                .collect(Collectors.toList());

        return limitBasedOnSupport(AssemblyFiltering.trimAndDeduplicate(mSupportChecker, rightExtended), 4);
    }

    private String nextAssemblyName()
    {
        return String.format("%sE%s", mPrimary.Name, mNextAssemblyNumber++);
    }

    private List<Pair<Record, Record>> findMates(final Collection<Record> records)
    {
        return mSAMSource.streamMates(records).collect(Collectors.toList());
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

    private HeadNode alignmentAsGraph(final Record alignment, final Direction orientation)
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

    private List<ExtendedAssembly> extendLeft(final SupportedAssembly assembly, final List<Record> discordantReads)
    {
        final List<Pair<Record, Record>> pairedMates = findMates(assembly.getSupportRecords().stream()
                .filter(Record::isMateOnTheLeft)
                .collect(Collectors.toList()));

        final List<Record> mates = pairedMates.stream()
                .filter(pair -> !assembly.containsSupport(pair.getRight()))
                .filter(pair -> AlignmentFilters.isRecordAverageQualityAbove(pair.getRight(), SvConstants.AVG_BASE_QUAL_THRESHOLD))
                .map(pair -> pair.getLeft().isPositiveStrand() == pair.getRight().isPositiveStrand()
                        ? pair.getRight().flipStrand()
                        : pair.getRight())
                .sorted(Comparator.comparingInt(Record::getUnclippedEnd).reversed())
                .collect(Collectors.toList());
        mCounters.LeftMates.add(mates.size());

        discordantReads.sort(Comparator.comparingInt(Record::getUnclippedEnd).reversed());
        final List<Record> applicableDiscordantReads = new ArrayList<>();
        for(Record record : discordantReads)
            if(!record.isMateOnTheLeft())
                applicableDiscordantReads.add(record);

        final Map<Record, Integer> checkStartIndices = new LinkedHashMap<>();
        for(Record record : Iterables.concat(mates, applicableDiscordantReads))
        {
            // Get the support index of this record or its mate
            @Nullable
            final Integer existingSupportIndex = assembly.getSupport(record.getName()).stream()
                    .map(entry -> entry.getValue() + entry.getKey().getLength())
                    .min(Comparator.naturalOrder())
                    .orElse(null);
            final int minDepth = existingSupportIndex == null
                    ? assembly.getLength() - record.getLength()
                    : Math.max(assembly.getLength() - record.getLength(), assembly.getLength() - existingSupportIndex);
            checkStartIndices.put(record, minDepth);
        }

        final Direction assemblyReadDirection = Direction.REVERSE;
        final Direction mateReadDirection = Direction.REVERSE;
        final List<ExtendedAssembly> extended = extendAssembly(assembly, assemblyReadDirection, mates, discordantReads,
                checkStartIndices, mateReadDirection, mCounters.LeftMatesAssembled, mCounters.LeftDiscordantReadsAssembled);
        mCounters.ExtendLeftAssemblies.add(extended.size());
        return extended;
    }

    private List<ExtendedAssembly> extendRight(final SupportedAssembly assembly, final List<Record> discordantReads)
    {
        final List<Pair<Record, Record>> pairedMates = findMates(assembly.getSupportRecords().stream()
                .filter(record -> !record.isMateOnTheLeft())
                .collect(Collectors.toList()));

        final List<Record> mates = pairedMates.stream()
                .filter(pair -> !assembly.containsSupport(pair.getRight()))
                .filter(pair -> AlignmentFilters.isRecordAverageQualityAbove(pair.getRight(), SvConstants.AVG_BASE_QUAL_THRESHOLD))
                .sorted(Comparator.comparingInt(pair -> assembly.getSupportIndex(pair.getLeft())))
                .map(pair -> pair.getLeft().isPositiveStrand() == pair.getRight().isPositiveStrand()
                        ? pair.getRight().flipStrand()
                        : pair.getRight())
                .collect(Collectors.toList());
        mCounters.RightMates.add(mates.size());

        discordantReads.sort(Comparator.comparingInt(Record::getUnclippedStart));
        final List<Record> applicableDiscordantReads = new ArrayList<>();
        for(Record record : discordantReads)
            if(record.isMateOnTheLeft())
                applicableDiscordantReads.add(record);

        final Map<Record, Integer> checkStartIndices = new LinkedHashMap<>();
        for(Record record : Iterables.concat(mates, applicableDiscordantReads))
        {
            // Get the support index of this record or its mate
            @Nullable
            final Integer existingSupportIndex = assembly.getSupport(record.getName()).stream()
                    .map(Map.Entry::getValue)
                    .max(Comparator.naturalOrder())
                    .orElse(null);
            final int minDepth = existingSupportIndex == null
                    ? assembly.getLength() - record.getLength()
                    : Math.max(assembly.getLength() - record.getLength(), existingSupportIndex);
            checkStartIndices.put(record, minDepth);
        }

        final Direction assemblyReadDirection = Direction.FORWARDS;
        final Direction mateReadDirection = Direction.FORWARDS;
        final List<ExtendedAssembly> extended = extendAssembly(assembly, assemblyReadDirection, mates, applicableDiscordantReads,
                checkStartIndices, mateReadDirection, mCounters.RightMatesAssembled, mCounters.RightDiscordantReadsAssembled);
        mCounters.ExtendRightAssemblies.add(extended.size());
        return extended;
    }

    private List<ExtendedAssembly> extendAssembly(final SupportedAssembly assembly, final Direction assemblyDirection,
            final List<Record> mateAlignments,
            final List<Record> discordantAlignments,
            final Map<Record, Integer> alignmentMinDepth,
            final Direction alignmentDirection,
            final Counter mateAttachCounter,
            final Counter discordantAttachCounter)
    {
        if(mateAlignments.isEmpty() && discordantAlignments.isEmpty())
        {
            if(assembly instanceof ExtendedAssembly)
                return List.of((ExtendedAssembly) assembly);
            final ExtendedAssembly newAssembly = new ExtendedAssembly(assembly.Name, assembly.Assembly, assembly);
            assembly.getSupport().forEach(entry -> newAssembly.addEvidenceAt(entry.getKey(), entry.getValue()));
            return List.of(newAssembly);
        }

        final HeadNode existing = HeadNode.create(assembly, assemblyDirection);

        final List<Record> potentialNewSupport = new ArrayList<>();
        extendAssembly(existing, potentialNewSupport, mateAlignments, alignmentMinDepth, alignmentDirection, mateAttachCounter);
        extendAssembly(existing, potentialNewSupport, discordantAlignments, alignmentMinDepth, alignmentDirection, discordantAttachCounter);

        final Map<Record, Set<Integer>> supportStartIndices = potentialSupportIndices(existing);

        final String diagramSetName = "Extension " + assemblyDirection.name().toLowerCase();
        @Nullable
        final DiagramSet diagrams = simplifyGraph(diagramSetName, existing);

        return existing.flatten().stream()
                .map(newAssembly ->
                {
                    final String directionCorrected = assemblyDirection == Direction.FORWARDS
                            ? newAssembly
                            : new StringBuilder(newAssembly).reverse().toString();
                    final ExtendedAssembly newCandidate = new ExtendedAssembly(nextAssemblyName(), directionCorrected, assembly);
                    newCandidate.addDiagrams(diagrams);

                    reAddSupport(newCandidate, assembly, potentialNewSupport, supportStartIndices,
                            assemblyDirection == Direction.FORWARDS);

                    return newCandidate;
                }).collect(Collectors.toList());
    }

    private void extendAssembly(
            final HeadNode existing, final List<Record> potentialNewSupport, final List<Record> alignments,
            final Map<Record, Integer> alignmentMinDepth, final Direction alignmentDirection, final Counter attachCounter)
    {
        for(final Record record : alignments)
        {
            @Nullable
            final Integer depth = alignmentMinDepth.get(record);
            final HeadNode toAttach = alignmentAsGraph(record, alignmentDirection);

            if(existing.attach(toAttach, SvConstants.ASSEMBLYEXTENSIONMINMATCHEDBASES, SvConstants.ASSEMBLYEXTENSIONMAXMISMATCHES,
                    Objects.requireNonNullElse(depth, 0), Integer.MAX_VALUE))
            {
                potentialNewSupport.add(record);
                attachCounter.add(1);
            }
        }
    }

    /** Once we've built the graph, try and work out where in the flattened output different records may appear. */
    private Map<Record, Set<Integer>> potentialSupportIndices(final HeadNode head)
    {
        final Map<Record, Set<Integer>> candidateLocations = new HashMap<>();
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
                    candidateLocations.computeIfAbsent(support.Record, r -> new LinkedHashSet<>()).add(depth - support.ReadIndex);

            for(Node successor : node.successors())
                worklist.add(Pair.of(successor, depth + 1));
        }

        return candidateLocations;
    }

    private void reAddSupport(final ExtendedAssembly assembly, final SupportedAssembly original, final List<Record> potentialSupport,
            final Map<Record, Set<Integer>> supportStartIndices, final boolean isForwards)
    {
        if(potentialSupport.size() > 1_000)
            return;

        final int firstIndex = assembly.Assembly.indexOf(original.Assembly);
        final int lastIndex = assembly.Assembly.indexOf(original.Assembly, firstIndex);
        if(firstIndex == -1 || firstIndex != lastIndex)
            addSupportNoOffset(assembly, original);
        else
            addSupportWithOffset(assembly, original, firstIndex);

        for(Record record : potentialSupport)
        {
            if(assembly.containsSupport(record))
                continue;

            @Nullable
            final Set<Integer> candidates = supportStartIndices.get(record);

            if(candidates == null)
            {
                assembly.tryAddSupport(mSupportChecker, record);
            }
            else
            {
                for(int rawIndex : candidates)
                {
                    final int checkIndex = isForwards
                            ? rawIndex
                            : assembly.getLength() - rawIndex - record.getLength();
                    if(mSupportChecker.WeakSupport.supportsAt(assembly, record, checkIndex))
                    {
                        assembly.addEvidenceAt(record, checkIndex);
                        break;
                    }
                }
            }
        }
    }

    private void addSupportWithOffset(final ExtendedAssembly assembly, final SupportedAssembly original, final int offset)
    {
        original.getSupport().forEach(entry ->
        {
            if(mSupportChecker.WeakSupport.supportsAt(assembly, entry.getKey(), entry.getValue() + offset))
                assembly.addEvidenceAt(entry.getKey(), entry.getValue() + offset);
            else
            {
                mSupportChecker.supportsAtIndex(original, entry.getKey(), 20, 2, entry.getValue());
                mSupportChecker.supportsAtIndex(original, entry.getKey(), 20, 2, -32);
            }
        });
    }

    private void addSupportNoOffset(final ExtendedAssembly assembly, final SupportedAssembly original)
    {
        original.getSupport().forEach(entry -> assembly.tryAddSupport(mSupportChecker, entry.getKey()));
    }
}
