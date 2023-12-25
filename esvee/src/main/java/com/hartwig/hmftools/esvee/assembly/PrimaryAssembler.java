package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.Direction;
import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.JunctionProcessingException;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.models.DiagramSet;
import com.hartwig.hmftools.esvee.models.PrimaryAssembly;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.sam.SAMSource;
import com.hartwig.hmftools.esvee.util.Counter;
import com.hartwig.hmftools.esvee.processor.PrimaryAssemblyResult;
import com.hartwig.hmftools.esvee.processor.Problem;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class PrimaryAssembler
{
    private final SvConfig mConfig;
    private final SAMSource mSAMSource;
    private final SupportChecker mSupportChecker;
    private final NodeFolder mNodeFolder;
    private final String mJunctionChromosome;
    private final int mJunctionPosition;
    private final Direction mJunctionOrientation;
    private final PrimaryAssemblerCounters mCounters = new PrimaryAssemblerCounters();
    private final boolean mCreateDiagrams;

    private int mNextAssemblyNumber = 1;

    @Nullable
    public static PrimaryAssemblyResult process(final Context context, final Junction junction, final PrimaryAssemblerCounters counters)
    {
        try
        {
            boolean createDiagrams = context.Config.writeHtmlFiles() && context.Config.PlotDiagrams;

            final PrimaryAssembler assembler = new PrimaryAssembler(context.Config, context.SAMSource,
                    context.SupportChecker, junction.chromosome(), junction.position(), junction.orientation(),
                    createDiagrams);

            final List<PrimaryAssembly> candidateAssemblies = assembler.processJunction(junction);
            counters.add(assembler.getCounters());

            return new PrimaryAssemblyResult(junction, counters, List.of(), candidateAssemblies, List.of());
        }
        catch(final Throwable throwable)
        {
            final Problem problem = new Problem("Failure during primary assembly", throwable, junction);
            context.Problems.add(problem);
            
            if(throwable instanceof JunctionProcessingException)
                SV_LOGGER.warn("{}", problem);
            else
                SV_LOGGER.error("{}", problem, throwable);
            return null;
        }
    }

    public PrimaryAssembler(final SvConfig config, final SAMSource source,
            final SupportChecker supportChecker, final String junctionChromosome, final int junctionPosition,
            final Direction orientation, final boolean createDiagrams)
    {
        mConfig = config;
        mSAMSource = source;
        mSupportChecker = supportChecker;
        mNodeFolder = new NodeFolder();
        mJunctionChromosome = junctionChromosome;
        mJunctionPosition = junctionPosition;
        mJunctionOrientation = orientation;
        mCreateDiagrams = createDiagrams;
    }

    public PrimaryAssemblerCounters getCounters()
    {
        return mCounters;
    }

    public List<PrimaryAssembly> processJunction(final Junction junction)
    {
        SV_LOGGER.trace("Processing {} junction @ {}:{}", junction.orientation().name().toLowerCase(),
                junction.chromosome(), junction.position());

        return doProcessJunction(junction);
    }

    private List<PrimaryAssembly> doProcessJunction(final Junction junction)
    {
        final List<Record> nearbyAlignments = mSAMSource.findReadsNear(junction.chromosome(), junction.position());

        final List<Record> rawAlignments = nearbyAlignments.stream()
                .filter(Counter.asPredicate(alignment -> AlignmentFilters.alignmentCrossesJunction(alignment, junction), mCounters.ReadsCrossingJunction))
                .filter(Counter.asPredicate(alignment -> AlignmentFilters.isRecordAverageQualityAbove(alignment, SvConstants.AVG_BASE_QUAL_THRESHOLD), mCounters.ReadsPassingRawQualityThreshold))
                .map(alignment -> realignForJunction(alignment, junction))
                .collect(Collectors.toList());

        final List<Record> withLowQAlignments = rawAlignments.stream()
                .filter(Counter.asPredicate(alignment -> AlignmentFilters.recordSoftClipsNearJunction(alignment, junction), mCounters.ReadsSoftClippedAtJunction))
                .collect(Collectors.toList());

        final List<Record> filteredAlignments = withLowQAlignments.stream()
                .filter(Counter.asPredicate(alignment -> AlignmentFilters.isRecordAverageQualityPastJunctionAbove(alignment, junction, SvConstants.AVG_BASE_QUAL_THRESHOLD), mCounters.ReadsPassingJunctionQualityThreshold))
                .filter(Counter.asPredicate(alignment -> AlignmentFilters.hasAcceptableMapQ(alignment, SvConstants.MINMAPQTOSTARTJUNCTION), mCounters.HasAcceptableMapQ))
                .filter(Counter.asPredicate(AlignmentFilters::isNotBadlyMapped, mCounters.WellMapped))
                .collect(Collectors.toList());

        if(filteredAlignments.isEmpty())
            return List.of(); // There are no reads of acceptable quality supporting this junction

        final List<PrimaryAssembly> initialAssemblies = createInitialAssemblies(filteredAlignments);

        mCounters.InitialAssemblies.add(initialAssemblies.size());

        final List<PrimaryAssembly> extendedInitial = SvConstants.EXTEND_PRIMARIES
                 ? extendInitial(withLowQAlignments, initialAssemblies) : initialAssemblies;

        final List<PrimaryAssembly> dedupedInitial = AssemblyFiltering.trimAndDeduplicate(mSupportChecker, extendedInitial);
        mCounters.DedupedInitialAssemblies.add(dedupedInitial.size());
        
        final List<PrimaryAssembly> anchored = createAnchors(rawAlignments, dedupedInitial);

        for(PrimaryAssembly assembly : anchored)
        {
            for(Record alignment : withLowQAlignments)
            {
                if(!assembly.containsSupport(alignment))
                {
                    @Nullable
                    final Integer supportIndex = mSupportChecker.WeakSupport.bestSupportIndex(assembly, alignment, 50);
                    if(supportIndex != null)
                        assembly.addEvidenceAt(alignment, supportIndex);
                }
            }
        }

        final List<PrimaryAssembly> assemblies = AssemblyFiltering.trimAndDeduplicate(mSupportChecker, anchored);
        mCounters.DedupedAnchoredAssemblies.add(assemblies.size());

        final JunctionMetrics junctionMetrics = new JunctionMetrics(mJunctionChromosome, mJunctionPosition, mJunctionOrientation, mCounters);
        assemblies.forEach(assembly -> assembly.addErrata(junctionMetrics));
        return assemblies;
    }

    private List<PrimaryAssembly> extendInitial(final List<Record> alignments, final List<PrimaryAssembly> assemblies)
    {
        return assemblies.stream()
                .map(assembly -> extendInitial(alignments, assembly, mJunctionOrientation))
                .collect(Collectors.toList());
    }

    /** There may be alignments that can extend the assembly but that are too noisy to be used during initial construction.
     * Examples of these types of alignments may be, for example, ones with larger soft-clips that have resulted in unacceptably low MapQ.
     * Extension in this manner is not supposed to create new candidates, so we will always choose the "best" result after pruning. */
    private PrimaryAssembly extendInitial(final List<Record> alignments, final PrimaryAssembly assembly, final Direction direction)
    {
        HeadNode graph = HeadNode.create(assembly, direction);
        final Set<Record> support = assembly.getSupportRecords().stream()
                .collect(Collectors.toSet());
        for(Record alignment : alignments)
        {
            if (support.contains(alignment))
                continue;
            if (!mSupportChecker.WeakSupport.supports(assembly, alignment))
                continue; // PERF: This should be supports-at

            graph = HeadNode.combine(graph, HeadNode.create(alignment, assembly.AnchorPosition, direction));
        }

        final var diagrams = simplifyGraph("Initial Extension", graph, true);
        final List<String> flattened = graph.flatten();
        if (flattened.size() * alignments.size() > 100_000)
            //throw new JunctionProcessingException("Too many flattened assemblies or alignments!");
            SV_LOGGER.info("{} got {} extensions & {} alignments for a product of {}",
                    assembly.getName(), flattened.size(), alignments.size(), flattened.size() * alignments.size());
        return Stream.concat(Stream.of(assembly), flattened.stream()
                .map(assemblyBases ->
                {
                    if (direction == Direction.REVERSE)
                        assemblyBases = new StringBuilder(assemblyBases).reverse().toString();

                    final int anchorPositionInAssembly = direction == Direction.FORWARDS
                            ? 1
                            : assemblyBases.length() - 1;
                    final PrimaryAssembly newAssembly = new PrimaryAssembly(nextAssemblyName(), assemblyBases, assembly.AnchorChromosome,
                                    assembly.AnchorPosition, anchorPositionInAssembly);
                    newAssembly.Diagrams.addAll(assembly.Diagrams);
                    newAssembly.Diagrams.add(diagrams);
                    for(Record record : alignments)
                    {
                        // To support the assembly we need to either be fully contained in the assembly, or to support
                        // it with our back half if we're a forwards junction / front half if we're a backwards junction.
                        final int minSupportIndex = mJunctionOrientation == Direction.FORWARDS
                                ? -record.getLength()
                                : 0;
                        final int maxSupportIndex = mJunctionOrientation == Direction.FORWARDS
                                ? Math.min(0, newAssembly.Assembly.length() - record.getLength())
                                : newAssembly.Assembly.length();

                        @Nullable
                        final Integer supportIndex = mSupportChecker.StrongSupport.supportIndex(newAssembly, record, 3, minSupportIndex, maxSupportIndex);
                        if(supportIndex != null)
                            newAssembly.addEvidenceAt(record, supportIndex);
                    }
                    return newAssembly;
                }))
                .max(Comparator.comparingInt(newAssembly -> newAssembly.getSupportFragments().size()))
                .orElseThrow();
    }

    private String nextAssemblyName()
    {
        return String.format("%s:%s%s:%s", mJunctionChromosome, mJunctionPosition,
                mJunctionOrientation == Direction.FORWARDS ? "F" : "R", mNextAssemblyNumber++);
    }

    /**
     * Convert indels near the junction to soft-clips
     */
    private Record realignForJunction(final Record alignment, final Junction junction)
    {
        boolean justHadIndel = false;
        boolean wasRightNearJunction = false;
        int referencePosition = alignment.getAlignmentStart();
        for(int i = 0; i < alignment.getCigar().numCigarElements(); i++)
        {
            final CigarElement element = alignment.getCigar().getCigarElement(i);

            final boolean isIndel = element.getOperator() == CigarOperator.D || element.getOperator() == CigarOperator.I;
            final boolean isLeftNearJunction = Math.abs(referencePosition - junction.position()) <= 2;
            final boolean isRightNearJunction = Math.abs(referencePosition + element.getLength() - junction.position()) <= 2;

            if(junction.orientation() == Direction.REVERSE && isLeftNearJunction && justHadIndel)
            {
                // Soft-clip everything before this cigar element (excl)
                int softClippedLength = 0;
                for(int j = 0; j < i; j++)
                    if(alignment.getCigar().getCigarElement(j).getOperator().consumesReadBases())
                        softClippedLength += alignment.getCigar().getCigarElement(j).getLength();

                final List<CigarElement> newElements = new ArrayList<>();
                newElements.add(new CigarElement(softClippedLength, CigarOperator.S));
                for(int j = i; j < alignment.getCigar().numCigarElements(); j++)
                    newElements.add(alignment.getCigar().getCigarElement(j));

                final Record newAlignment = alignment.copy();
                newAlignment.setAlignmentStart(referencePosition);
                newAlignment.setCigar(new Cigar(newElements));
                return newAlignment;
            }

            if(junction.orientation() == Direction.FORWARDS && wasRightNearJunction && isIndel)
            {
                // Soft-clip everything after this cigar element (incl)
                int softClippedLength = 0;
                for(int j = i; j < alignment.getCigar().numCigarElements(); j++)
                    if(alignment.getCigar().getCigarElement(j).getOperator().consumesReadBases())
                        softClippedLength += alignment.getCigar().getCigarElement(j).getLength();

                final List<CigarElement> newElements = new ArrayList<>();
                for(int j = 0; j < i; j++)
                    newElements.add(alignment.getCigar().getCigarElement(j));
                newElements.add(new CigarElement(softClippedLength, CigarOperator.S));

                final Record newAlignment = alignment.copy();
                newAlignment.setCigar(new Cigar(newElements));
                return newAlignment;
            }

            if(element.getOperator().consumesReferenceBases())
                referencePosition += element.getLength();
            justHadIndel = isIndel;
            wasRightNearJunction = isRightNearJunction;
        }

        return alignment;
    }

    private List<PrimaryAssembly> createInitialAssemblies(final List<Record> alignments)
    {
        final HeadNode combinedForwards = alignments.stream()
                .filter(alignment -> alignment.getChromosome().equals(mJunctionChromosome))
                .map(alignment -> HeadNode.create(alignment, mJunctionPosition, mJunctionOrientation))
                .filter(Objects::nonNull)
                .reduce(HeadNode::combine)
                .orElseThrow();

        @Nullable
        final DiagramSet diagrams = simplifyGraph("Initial Construction", combinedForwards, false);

        final List<String> flattened = combinedForwards.flatten();

        mCounters.FlattenedInitial.add(flattened.size());
        final List<PrimaryAssembly> candidateAssemblies = flattened.stream()
                .filter(assemblyString -> assemblyString.length() >= 10)
                .map(assemblyString ->
                {
                    final String orientedAssembly = mJunctionOrientation == Direction.FORWARDS
                            ? assemblyString
                            : new StringBuilder(assemblyString).reverse().toString();
                    final int anchorPositionInAssembly = mJunctionOrientation == Direction.FORWARDS
                            ? 0
                            : orientedAssembly.length() - 1;

                    return new PrimaryAssembly(nextAssemblyName(), orientedAssembly, mJunctionChromosome, mJunctionPosition, anchorPositionInAssembly);
                })
                .peek(assembly -> assembly.addDiagrams(diagrams))
                .collect(Collectors.toList());

        for(PrimaryAssembly assembly : candidateAssemblies)
        {
            for(Record record : alignments)
            {
                // To support the assembly we need to either be fully contained in the assembly, or to support
                // it with our back half if we're a forwards junction / front half if we're a backwards junction.
                final int minSupportIndex = mJunctionOrientation == Direction.FORWARDS
                        ? -record.getLength()
                        : 0;
                final int maxSupportIndex = mJunctionOrientation == Direction.FORWARDS
                        ? Math.min(0, assembly.Assembly.length() - record.getLength())
                        : assembly.Assembly.length();

                @Nullable
                final Integer supportIndex = mSupportChecker.WeakSupport.supportIndex(assembly, record, 3, minSupportIndex, maxSupportIndex);
                if(supportIndex != null)
                    assembly.addEvidenceAt(record, supportIndex);
            }
        }

        return candidateAssemblies.stream()
                .filter(assembly -> assembly.supportCount() != 0)
                .collect(Collectors.toList());
    }

    private DiagramSet simplifyGraph(final String diagramSetName, final HeadNode node, final boolean aggressive)
    {
        final DiagramSet diagrams = new DiagramSet(diagramSetName);

        if(mCreateDiagrams)
            diagrams.add("Attachment", node.toDiagram());

        mNodeFolder.foldPaths(node);

        if(mCreateDiagrams)
            diagrams.add("Folding", node.toDiagram());

        if(aggressive)
            node.pruneNodesAggressive();
        else
            node.pruneNodes();

        if(mCreateDiagrams)
            diagrams.add("Pruning", node.toDiagram());

        return diagrams;
    }

    private List<PrimaryAssembly> createAnchors(final List<Record> alignments, final List<PrimaryAssembly> initialAssemblies)
    {
        final Map<Record, HeadNode> reverseSequences = new HashMap<>();
        for(Record record : alignments)
            reverseSequences.put(record, HeadNode.create(record, mJunctionPosition, mJunctionOrientation.opposite()));

        final List<PrimaryAssembly> anchored = initialAssemblies.stream()
                .flatMap(candidateAssembly -> createAnchor(reverseSequences, candidateAssembly).stream())
                .collect(Collectors.toList());

        mCounters.AnchoredAssemblies.add(anchored.size());
        return anchored;
    }

    private List<PrimaryAssembly> createAnchor(final Map<Record, HeadNode> reverseSequences, final PrimaryAssembly initialAssembly)
    {
        HeadNode anchor = null;
        for(Record support : initialAssembly.getSupportRecords())
        {
            final HeadNode sequence = reverseSequences.get(support);
            if(anchor == null)
                anchor = sequence.deepCopy();
            else
                anchor = HeadNode.combine(anchor, sequence, false);
        }
        assert anchor != null;
        anchor.sortSupport();
        anchor = anchor.deepCopy();

        @Nullable
        final DiagramSet diagrams = simplifyGraph("Anchor Construction", anchor, false);

        final List<String> flattened = anchor.flatten();
        mCounters.FlattenedAnchors.add(flattened.size());
        final List<PrimaryAssembly> anchoredAssemblies = new ArrayList<>();
        for(String flattenedAssembly : flattened)
        {
            final boolean isForwards = mJunctionOrientation == Direction.FORWARDS;
            final String anchoredAssembly = isForwards
                    ? new StringBuilder(flattenedAssembly.substring(1)).reverse() + initialAssembly.Assembly
                    : initialAssembly.Assembly + flattenedAssembly.substring(1);

            final int offsetSize = isForwards ? flattenedAssembly.length() - 1 : 0;
            final int anchorPositionInAssembly = initialAssembly.AnchorPositionInAssembly + offsetSize;

            // Don't use candidate to construct, as we want to re-evaluate evidence
            final PrimaryAssembly assembly = new PrimaryAssembly(nextAssemblyName(), anchoredAssembly,
                    mJunctionChromosome, mJunctionPosition, anchorPositionInAssembly);
            assembly.Diagrams.addAll(initialAssembly.Diagrams);
            assembly.addDiagrams(diagrams);

            for(var entry : initialAssembly.getSupport())
            {
                final Record record = entry.getKey();
                final int supportIndex = entry.getValue();
                final int newSupportIndex;
                if (isForwards)
                    newSupportIndex = supportIndex + (anchoredAssembly.length() - initialAssembly.Assembly.length());
                else
                    newSupportIndex = supportIndex;

                if (mSupportChecker.WeakSupport.supportsAt(assembly, record, newSupportIndex))
                    assembly.addEvidenceAt(record, newSupportIndex);
                else
                    assembly.tryAddSupport(mSupportChecker, record);
            }

            if (assembly.getSupportFragments().size() > SvConstants.MINREADSTOSUPPORTASSEMBLY)
                anchoredAssemblies.add(assembly);
        }

        return anchoredAssemblies;
    }
}
