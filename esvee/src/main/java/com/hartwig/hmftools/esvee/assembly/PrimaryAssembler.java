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

import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.html.DiagramSet;
import com.hartwig.hmftools.esvee.sequence.PrimaryAssembly;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.util.Counter;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class PrimaryAssembler
{
    private final SvConfig mConfig;

    private final Junction mJunction;
    
    private final SupportChecker mSupportChecker;
    private final NodeFolder mNodeFolder;
    private final PrimaryAssemblerCounters mCounters;
    private final boolean mCreateDiagrams;

    private int mNextAssemblyNumber = 1;

    public PrimaryAssembler(final SvConfig config, final Junction junction)
    {
        mConfig = config;
        mSupportChecker = new SupportChecker();
        mNodeFolder = new NodeFolder();
        mCounters = new PrimaryAssemblerCounters();
        mJunction = junction;
        mCreateDiagrams = config.writeHtmlFiles() && config.PlotDiagrams;
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
                .filter(Counter.asPredicate(alignment -> AlignmentFilters.recordSoftClipsNearJunction(alignment, mJunction), mCounters.ReadsSoftClippedAtJunction))
                .collect(Collectors.toList());

        final List<Read> filteredAlignments = withLowQAlignments.stream()
                .filter(Counter.asPredicate(alignment -> AlignmentFilters.isRecordAverageQualityPastJunctionAbove(alignment, mJunction, SvConstants.AVG_BASE_QUAL_THRESHOLD), mCounters.ReadsPassingJunctionQualityThreshold))
                .filter(Counter.asPredicate(alignment -> AlignmentFilters.hasAcceptableMapQ(alignment, SvConstants.MIN_MAPQ_START_JUNCTION), mCounters.HasAcceptableMapQ))
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
        
        final List<PrimaryAssembly> anchored = createAnchors(realignedReads, dedupedInitial);

        for(PrimaryAssembly assembly : anchored)
        {
            for(Read alignment : withLowQAlignments)
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

        final JunctionMetrics junctionMetrics = new JunctionMetrics(mJunction.Chromosome, mJunction.Position, mJunction.direction(), mCounters);
        assemblies.forEach(assembly -> assembly.addErrata(junctionMetrics));
        return assemblies;
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
    private PrimaryAssembly extendInitial(final List<Read> alignments, final PrimaryAssembly assembly, final Direction direction)
    {
        HeadNode graph = HeadNode.create(assembly, direction);
        final Set<Read> support = assembly.getSupportRecords().stream()
                .collect(Collectors.toSet());
        for(Read alignment : alignments)
        {
            if(support.contains(alignment))
                continue;
            if(!mSupportChecker.WeakSupport.supports(assembly, alignment))
                continue; // PERF: This should be supports-at

            graph = HeadNode.combine(graph, HeadNode.create(alignment, assembly.AnchorPosition, direction));
        }

        final var diagrams = simplifyGraph("Initial Extension", graph, true);
        final List<String> flattened = graph.flatten();
        if(flattened.size() * alignments.size() > 100_000)
            //throw new JunctionProcessingException("Too many flattened assemblies or alignments!");
            SV_LOGGER.info("{} got {} extensions & {} alignments for a product of {}",
                    assembly.getName(), flattened.size(), alignments.size(), flattened.size() * alignments.size());
        return Stream.concat(Stream.of(assembly), flattened.stream()
                .map(assemblyBases ->
                {
                    if(direction == Direction.REVERSE)
                        assemblyBases = new StringBuilder(assemblyBases).reverse().toString();

                    final int anchorPositionInAssembly = direction == Direction.FORWARDS
                            ? 1
                            : assemblyBases.length() - 1;
                    final PrimaryAssembly newAssembly = new PrimaryAssembly(nextAssemblyName(), assemblyBases, assembly.AnchorChromosome,
                                    assembly.AnchorPosition, anchorPositionInAssembly);
                    newAssembly.Diagrams.addAll(assembly.Diagrams);
                    newAssembly.Diagrams.add(diagrams);
                    for(Read read : alignments)
                    {
                        // To support the assembly we need to either be fully contained in the assembly, or to support
                        // it with our back half if we're a forwards junction / front half if we're a backwards junction.
                        final int minSupportIndex = mJunction.direction() == Direction.FORWARDS
                                ? -read.getLength()
                                : 0;
                        final int maxSupportIndex = mJunction.direction() == Direction.FORWARDS
                                ? Math.min(0, newAssembly.Assembly.length() - read.getLength())
                                : newAssembly.Assembly.length();

                        @Nullable
                        final Integer supportIndex = mSupportChecker.StrongSupport.supportIndex(newAssembly, read, 3, minSupportIndex, maxSupportIndex);
                        if(supportIndex != null)
                            newAssembly.addEvidenceAt(read, supportIndex);
                    }
                    return newAssembly;
                }))
                .max(Comparator.comparingInt(newAssembly -> newAssembly.getSupportFragments().size()))
                .orElseThrow();
    }

    private String nextAssemblyName()
    {
        return String.format("%s:%s%s:%s", mJunction.Chromosome, mJunction.Position,
                mJunction.direction() == Direction.FORWARDS ? "F" : "R", mNextAssemblyNumber++);
    }

    // convert indels near the junction to soft-clips
    private Read realignForJunction(final Read read, final Junction junction)
    {
        boolean justHadIndel = false;
        boolean wasRightNearJunction = false;
        int referencePosition = read.getAlignmentStart();
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

                // CHECK: decide how to handle these
                /*
                final Read newAlignment = read.copyRecord();
                newAlignment.setAlignmentStart(referencePosition);
                newAlignment.setCigar(new Cigar(newElements));
                return newAlignment;
                */
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
                .filter(alignment -> alignment.getChromosome().equals(mJunction.Chromosome))
                .map(alignment -> HeadNode.create(alignment, mJunction.Position, mJunction.direction()))
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
                    final String orientedAssembly = mJunction.direction() == Direction.FORWARDS
                            ? assemblyString
                            : new StringBuilder(assemblyString).reverse().toString();
                    final int anchorPositionInAssembly = mJunction.direction() == Direction.FORWARDS
                            ? 0
                            : orientedAssembly.length() - 1;

                    return new PrimaryAssembly(nextAssemblyName(), orientedAssembly, mJunction.Chromosome, mJunction.Position, anchorPositionInAssembly);
                })
                .peek(assembly -> assembly.addDiagrams(diagrams))
                .collect(Collectors.toList());

        for(PrimaryAssembly assembly : candidateAssemblies)
        {
            for(Read read : alignments)
            {
                // To support the assembly we need to either be fully contained in the assembly, or to support
                // it with our back half if we're a forwards junction / front half if we're a backwards junction.
                final int minSupportIndex = mJunction.direction() == Direction.FORWARDS
                        ? -read.getLength()
                        : 0;
                final int maxSupportIndex = mJunction.direction() == Direction.FORWARDS
                        ? Math.min(0, assembly.Assembly.length() - read.getLength())
                        : assembly.Assembly.length();

                @Nullable
                final Integer supportIndex = mSupportChecker.WeakSupport.supportIndex(assembly, read, 3, minSupportIndex, maxSupportIndex);
                if(supportIndex != null)
                    assembly.addEvidenceAt(read, supportIndex);
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

    private List<PrimaryAssembly> createAnchors(final List<Read> alignments, final List<PrimaryAssembly> initialAssemblies)
    {
        final Map<Read, HeadNode> reverseSequences = new HashMap<>();
        for(Read read : alignments)
            reverseSequences.put(read, HeadNode.create(read, mJunction.Position, mJunction.direction().opposite()));

        final List<PrimaryAssembly> anchored = initialAssemblies.stream()
                .flatMap(candidateAssembly -> createAnchor(reverseSequences, candidateAssembly).stream())
                .collect(Collectors.toList());

        mCounters.AnchoredAssemblies.add(anchored.size());
        return anchored;
    }

    private List<PrimaryAssembly> createAnchor(final Map<Read, HeadNode> reverseSequences, final PrimaryAssembly initialAssembly)
    {
        HeadNode anchor = null;
        for(Read support : initialAssembly.getSupportRecords())
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
            final boolean isForwards = mJunction.direction() == Direction.FORWARDS;
            final String anchoredAssembly = isForwards
                    ? new StringBuilder(flattenedAssembly.substring(1)).reverse() + initialAssembly.Assembly
                    : initialAssembly.Assembly + flattenedAssembly.substring(1);

            final int offsetSize = isForwards ? flattenedAssembly.length() - 1 : 0;
            final int anchorPositionInAssembly = initialAssembly.AnchorPositionInAssembly + offsetSize;

            // Don't use candidate to construct, as we want to re-evaluate evidence
            final PrimaryAssembly assembly = new PrimaryAssembly(nextAssemblyName(), anchoredAssembly,
                    mJunction.Chromosome, mJunction.Position, anchorPositionInAssembly);
            assembly.Diagrams.addAll(initialAssembly.Diagrams);
            assembly.addDiagrams(diagrams);

            for(var entry : initialAssembly.getSupport())
            {
                final Read read = entry.getKey();
                final int supportIndex = entry.getValue();
                final int newSupportIndex;
                if(isForwards)
                    newSupportIndex = supportIndex + (anchoredAssembly.length() - initialAssembly.Assembly.length());
                else
                    newSupportIndex = supportIndex;

                if(mSupportChecker.WeakSupport.supportsAt(assembly, read, newSupportIndex))
                    assembly.addEvidenceAt(read, newSupportIndex);
                else
                    assembly.tryAddSupport(mSupportChecker, read);
            }

            if(assembly.getSupportFragments().size() > SvConstants.MIN_READS_SUPPORT_ASSEMBLY)
                anchoredAssemblies.add(assembly);
        }

        return anchoredAssemblies;
    }
}
