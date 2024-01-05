package com.hartwig.hmftools.esvee.assembly;

import static htsjdk.samtools.CigarOperator.M;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.sequence.ReadSupport;
import com.hartwig.hmftools.esvee.sequence.SupportedAssembly;
import com.hartwig.hmftools.esvee.read.ReadUtils;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;

public class HeadNode extends Node
{
    private int mJunctions;

    public HeadNode()
    {
        super('S');
    }

    // approximate junction count, inaccurate after folding or mutating
    public int junctionCount()
    {
        return mJunctions;
    }

    public void trimShortPaths(final int threshold)
    {
        final Map<Node, Integer> depthLookup = new HashMap<>();
        final Deque<Pair<Node, Integer>> toProcess = new ArrayDeque<>();
        toProcess.add(Pair.of(this, -1));

        while(!toProcess.isEmpty())
        {
            final Pair<Node, Integer> pair = toProcess.removeLast();
            final Node node = pair.getLeft();
            final int length = pair.getRight();

            final Collection<Node> successors = node.successors();
            int nodeDepth = 1;
            for(Node successor : successors)
            {
                final int nextDepth = depthLookup.getOrDefault(successor, -1);
                nodeDepth = nextDepth == -1 ? -1 : Math.max(nextDepth + 1, nodeDepth);
                if(nodeDepth == -1)
                    break;

                if(nextDepth + length + 1 < threshold)
                    node.removeSuccessor(successor.Base);
            }

            if(nodeDepth == -1)
            {
                // Haven't computed successor yet
                toProcess.addLast(pair);
                successors.forEach(successor -> toProcess.addLast(Pair.of(successor, length + 1)));
            }
            else
                depthLookup.put(node, nodeDepth);
        }
    }

    private boolean isLowQualityBase(final Node node, final boolean checkDepth)
    {
        return node.MaxQuality <= SvConstants.LOW_BASE_QUAL_THRESHOLD
            || (checkDepth && node.supportDepth() > 1 && node.Quality <= SvConstants.LOW_BASE_QUAL_CUMULATIVE_THRESHOLD);
    }

    static List<Node.Support> mergeSupport(final List<Node.Support> left, final List<Node.Support> right)
    {
        if(left.size() == 1 && right.size() == 1 && left.get(0).compareTo(right.get(0)) == 0)
            return left;

        final List<Node.Support> merged = new ArrayList<>(left.size() + right.size());
        int leftIndex = 0;
        int rightIndex = 0;
        while(leftIndex < left.size() && rightIndex < right.size())
        {
            final Support nextLeft = left.get(leftIndex);
            final Support nextRight = right.get(rightIndex);

            final int comparison = nextLeft.compareTo(nextRight);
            if(comparison < 0)
            {
                merged.add(nextLeft);
                leftIndex++;
            }
            else if(comparison > 0)
            {
                merged.add(nextRight);
                rightIndex++;
            }
            else
            {
                merged.add(nextLeft);
                leftIndex++;
                rightIndex++;
            }
        }
        while(leftIndex < left.size())
            merged.add(left.get(leftIndex++));
        while(rightIndex < right.size())
            merged.add(right.get(rightIndex++));

        return merged;
    }

    public void pruneNodes()
    {
        pruneNodes(false);
    }

    public void pruneNodesAggressive()
    {
        pruneNodes(true);
    }

    private void pruneNodes(final boolean heavyHanded)
    {
        final Queue<Node> toProcess = new ArrayDeque<>();
        final Set<Node> visited = Collections.newSetFromMap(new IdentityHashMap<>());
        toProcess.add(this);

        while(!toProcess.isEmpty())
        {
            final Node node = toProcess.poll();
            if(!visited.add(node))
                continue;

            final var successors = node.successors();
            if(successors.size() == 1)
            {
                toProcess.add(successors.get(0));
                continue;
            }
            else if(successors.size() == 0)
                continue;

            // We have >1 successor
            // If any of our successors have been merged, determine if there's enough support for both bases
            for(int i = 0; i < successors.size(); i++)
                for(int j = 0; j < successors.size(); j++)
                {
                    final Node left = successors.get(i);
                    final Node right = successors.get(j);

                    if(i == j || left.Quality < right.Quality)
                        continue; // Look at from other side
                    assert left.Quality >= right.Quality;

                    final List<Node> leftSuccessors = left.successors();
                    final List<Node> rightSuccessors = right.successors();
                    //noinspection SlowListContainsAll
                    if(!leftSuccessors.containsAll(rightSuccessors))
                    {
                        if(rightSuccessors.isEmpty()
                                && right.Quality < SvConstants.LOW_BASE_QUAL_CUMULATIVE_THRESHOLD
                                && left.Quality > SvConstants.LOW_BASE_QUAL_CUMULATIVE_THRESHOLD)
                            node.removeSuccessor(right.Base);
                        continue; // Not merged
                    }
                    if(heavyHanded && right.Quality < left.Quality)
                    {
                        node.removeSuccessor(right.Base);
                        continue;
                    }

                    if(left.MaxQuality < SvConstants.LOW_BASE_QUAL_THRESHOLD && right.MaxQuality < SvConstants.LOW_BASE_QUAL_THRESHOLD)
                    {
                        if(left.Base < right.Base)
                        {
                            node.removeSuccessor(right.Base);
                            left.MaxQuality = left.Quality = 0;
                        }
                    }
                    else if(right.MaxQuality < SvConstants.LOW_BASE_QUAL_THRESHOLD)
                        node.removeSuccessor(right.Base);
                    else if(right.supportDepth() == 1 && left.supportDepth() > 2)
                        node.removeSuccessor(right.Base);
                    else if(right.supportDepth() == 2 && left.supportDepth() > 20)
                        node.removeSuccessor(right.Base);
                    else if(right.Quality == 0 && left.Quality != 0)
                        node.removeSuccessor(right.Base);
                }
            toProcess.addAll(node.successors());
        }
    }

    private static class FlattenEntry
    {
        public final Node Node;
        public final StringBuilder StringBuilder;
        public final int TotalQuality;

        private FlattenEntry(final Node node, final StringBuilder stringBuilder, final int totalQuality)
        {
            Node = node;
            StringBuilder = stringBuilder;
            TotalQuality = totalQuality;
        }

        public double averageQuality()
        {
            return (double) TotalQuality / StringBuilder.length();
        }
    }

    public List<String> flatten()
    {
        final List<String> results = new ArrayList<>();
        final Queue<FlattenEntry> toProcess = new ArrayDeque<>();
        toProcess.add(new FlattenEntry(this, new StringBuilder(), 0));
        while(!toProcess.isEmpty())
        {
            if(toProcess.size() > 1_000)
            {
                // Keep only the 900 entries with the "best" average quality
                final List<FlattenEntry> entries = new ArrayList<>(toProcess);
                entries.sort(Comparator.comparingDouble(FlattenEntry::averageQuality).reversed());

                toProcess.clear();
                for(int i = 0; i < 900; i++)
                    toProcess.add(entries.get(i));
                assert !toProcess.isEmpty();
            }

            final FlattenEntry entry = toProcess.poll();
            final Node currentNode = entry.Node;
            final StringBuilder sb = entry.StringBuilder;
            int totalQuality = entry.TotalQuality;

            if(!(currentNode instanceof HeadNode))
            {
                sb.append(currentNode.Base);
                totalQuality += currentNode.Quality;
            }

            final int successorCount = currentNode.successorCount();
            if(successorCount == 0 && sb.length() > 0)
                results.add(sb.toString());

            for(Node successor : currentNode.successors())
                if(successorCount == 1)
                    toProcess.add(new FlattenEntry(successor, sb, totalQuality));
                else
                    toProcess.add(new FlattenEntry(successor, new StringBuilder().append(sb), totalQuality));
        }

        return results;
    }

    private static class OverlayHead
    {
        public final Node Node;
        public final int Overlap;
        public final int Errors;
        public final Deque<Node> Path;

        private OverlayHead(final Node node, final int overlap, final int errors, final Deque<Node> path)
        {
            Node = node;
            Overlap = overlap;
            Errors = errors;
            Path = path;
        }

        private OverlayHead extend(final Node extension, final boolean isMatch, final boolean isMismatch, final boolean mutateHead)
        {
            final Deque<Node> extendedPath;
            if(mutateHead)
                extendedPath = Path;
            else
            {
                extendedPath = new ArrayDeque<>(Path.size() + 16);
                extendedPath.addAll(Path);
            }
            extendedPath.add(extension);
            return new OverlayHead(extension, Overlap + (isMatch ? 1 : 0), Errors + (isMismatch ? 1 : 0), extendedPath);
        }
    }

    @Nullable
    private Pair<Deque<Node>, Integer> scoreOverlay(final Node destination, final Node source, final int minimumOverlap, final int maxErrors)
    {
        // Do a breadth-first-search from destination, using source to rate how we do
        final Map<Node, OverlayHead> frozenHeads = new IdentityHashMap<>(4);
        Map<Node, OverlayHead> heads = new IdentityHashMap<>(4);
        final boolean isStartLowQuality = isLowQualityBase(destination, true) || isLowQualityBase(source, false);
        final boolean isStartMismatch = !isStartLowQuality && destination.Base != source.Base;

        final Deque<Node> newQueue = new ArrayDeque<>(16);
        newQueue.add(destination);
        heads.put(destination, new OverlayHead(destination, isStartMismatch ? 0 : 1, isStartMismatch ? 1 : 0, newQueue));

        // Sentinel for a path that is ambiguous
        final OverlayHead ambiguous = new OverlayHead(null, -1, 0, null);

        Node current = source;
        int bestOverlap = 0;
        while(!heads.isEmpty())
        {
            if(current.successors().isEmpty())
                break;
            current = current.successors().get(0);

            final Map<Node, OverlayHead> newHeads = new IdentityHashMap<>(4);
            for(OverlayHead head : heads.values())
            {
                if(head == ambiguous)
                    continue;

                final List<Node> successors = head.Node.successors();
                if(successors.isEmpty())
                {
                    if(head.Overlap >= bestOverlap)
                        frozenHeads.put(head.Node, head);
                    continue;
                }

                final boolean mutateHead = successors.size() == 1;
                for(Node successor : successors)
                {
                    final boolean isLowQuality = isLowQualityBase(successor, true) || isLowQualityBase(current, false);
                    final boolean isMatch = successor.Base == current.Base;
                    final boolean isMismatch = !isLowQuality && successor.Base != current.Base;
                    final int overlap = head.Overlap + (isMatch ? 1 : 0);
                    bestOverlap = Math.max(overlap, bestOverlap);

                    final int errors = head.Errors + (isMismatch ? 1 : 0);
                    if(errors > maxErrors)
                        continue;

                    final int currentQuality = current.MaxQuality;
                    newHeads.compute(successor, (n, existing) -> {
                        if(existing == null)
                            return head.extend(n, isMatch, isMismatch, mutateHead);
                        else if(existing == ambiguous)
                            return ambiguous;
                        else
                        {
                            if(existing.Overlap > overlap || (existing.Overlap == overlap && existing.Errors < errors))
                                return existing;
                            else if(existing.Overlap == overlap && existing.Errors == errors)
                            {
                                final Node otherCurrent = existing.Path.peekLast();
                                assert otherCurrent != null;
                                if(otherCurrent.MaxQuality == currentQuality)
                                    return ambiguous;
                                else if(otherCurrent.MaxQuality > currentQuality)
                                    return existing;
                                else
                                    return head.extend(n, isMatch, isMismatch, mutateHead);
                            }
                            else
                                return head.extend(n, isMatch, isMismatch, mutateHead);
                        }
                    });
                }
            }

            heads = newHeads;
        }
        heads.putAll(frozenHeads);

        if(heads.isEmpty())
            return null;

        final int bestScore = heads.values().stream()
                .mapToInt(head -> head.Overlap)
                .max().orElse(-1);
        if(bestScore < minimumOverlap)
            return null; // Does not meet minimum

        final List<OverlayHead> bestHeads = heads.values().stream()
                .filter(head -> head.Overlap == bestScore)
                .collect(Collectors.toList());
        if(bestHeads.size() != 1)
            return null; // Ambiguous
        final OverlayHead bestHead = bestHeads.get(0);
        final int neutralBases = bestHead.Path.size() - bestHead.Overlap;
        if(neutralBases >= bestHead.Overlap)
            return null; // Too many neutral bases

        return Pair.of(bestHead.Path, bestHead.Overlap);
    }

    private void doAttach(Node previous, final Deque<Node> destinationOverlay, final Node source)
    {
        @Nullable Node nextSource = source;
        while(!destinationOverlay.isEmpty())
        {
            final Node destination = destinationOverlay.poll();
            if(nextSource == null)
                throw new IllegalStateException("More to overlay, but no more source data?");

            if(nextSource.Base != destination.Base)
            {
                final Node node = new Node(nextSource.Base);
                node.MaxQuality = nextSource.MaxQuality;
                node.Quality = nextSource.Quality;
                node.Support = nextSource.Support;
                if(nextSource.successors().isEmpty())
                {
                    previous.setNext(node);
                    destination.successors().forEach(node::setNext);
                    return;
                }

                // A - T -  A
                //        \ T
                // +
                // A - G - A
                // ->
                //        / T
                // A - T -- A
                //   \   /
                //     G

                // A - T -  A
                //        \ T
                // +
                // A - G - G
                // ->
                //        / T
                // A - T -- A
                //   \
                //     G -- G
                //noinspection AssertWithSideEffects
                assert nextSource.successors().size() <= 1 : "Source had >1 successor";
                final Node singleSuccessor = nextSource.successors().get(0);
                final char nextBase = singleSuccessor.Base;
                @Nullable
                final Node destinationSuccessor = destination.successors().stream()
                        .filter(successor -> successor.Base == nextBase)
                        .findFirst().orElse(null);
                if(destinationSuccessor != null)
                {
                    if(node.MaxQuality > 12 || destinationSuccessor.MaxQuality < 35)
                    {
                        previous.setNext(node);
                        node.setNext(destinationSuccessor);
                    }
                }
                else
                {
                    // Linear copy from here
                    previous.setNext(node);
                    node.setNext(singleSuccessor.deepCopy());
                    return;
                }
            }
            else
            {
                destination.Support = mergeSupport(destination.Support, nextSource.Support);
                destination.recomputeQuality();
            }

            final List<Node> nextSuccessors = nextSource.successors();
            if(!nextSuccessors.isEmpty())
                nextSource = nextSuccessors.get(0);
            else
                nextSource = null;

            previous = destination;
        }

        if(nextSource != null)
            previous.setNext(nextSource);
    }

    public boolean attach(final HeadNode other, final int minimumOverlap, final int maxErrors, final int minDepth, final int maxDepth)
    {
        // locate the best place to overlay other onto this, and does so
        final List<Node> otherSuccessors = other.successors();
        if(otherSuccessors.size() != 1)
            throw new IllegalArgumentException("Can only attach simple sequences");
        final Node otherFirst = otherSuccessors.get(0);

        // At every position, score the overlay.
        int bestChildOverlap = -1;
        @Nullable
        Deque<Node> bestOverlay = null;
        Node previousForBestOverlay = this;

        // Triple<Node, Parent, DepthIndex>
        final Deque<Triple<Node, Node, Integer>> worklist = new ArrayDeque<>();
        successors().forEach(successor -> worklist.add(Triple.of(successor, this, 0)));

        // Required to deal with converging paths being identical
        final Set<Node> processed = Collections.newSetFromMap(new IdentityHashMap<>());
        while(!worklist.isEmpty())
        {
            final Triple<Node, Node, Integer> pair = worklist.poll();
            final Node node = pair.getLeft();
            final Node previous = pair.getMiddle();
            final int depth = pair.getRight();
            if(!processed.add(node) || depth > maxDepth)
                continue;

            if(depth >= minDepth)
            {
                @Nullable
                final Pair<Deque<Node>, Integer> overlayPair = scoreOverlay(node, otherFirst, minimumOverlap, maxErrors);
                @Nullable
                final Deque<Node> overlay = overlayPair == null ? null : overlayPair.getLeft();
                int overlap = overlayPair == null ? -1 : overlayPair.getRight();
                if(overlay != null && overlap >= bestChildOverlap)
                {
                    overlap = rescoreOverlay(overlay, overlap);

                    if(overlap > bestChildOverlap && overlap > minimumOverlap)
                    {
                        bestChildOverlap = overlap;
                        bestOverlay = overlay;
                        previousForBestOverlay = previous;
                    }
                    else if(overlap > 0 && overlap == bestChildOverlap)
                        return false; // Ambiguous overlay
                }
            }

            for(Node successor : node.successors())
            {
                worklist.add(Triple.of(successor, node, depth + 1));
            }
        }
        if(bestChildOverlap < minimumOverlap)
            return false;

        assert bestOverlay != null;
        doAttach(previousForBestOverlay, bestOverlay, otherFirst);

        return true;
    }

    private int rescoreOverlay(final Deque<Node> path, final int rawScore)
    {
        final int repeatLength = SvConstants.ASSEMBLY_EXTENSION_MAX_REPEAT_SCORE + 1;
        final int errors = path.size() - rawScore;
        int score = 0;
        final Deque<Node> recentPath = new ArrayDeque<>();
        for(Node n : path)
        {
            recentPath.addLast(n);
            if(recentPath.size() > repeatLength)
                recentPath.removeFirst();

            if(recentPath.size() != repeatLength || !recentPath.stream().allMatch(recent -> recent.Base == n.Base))
                score++;
        }

        return score - errors;
    }

    @Override
    public HeadNode deepCopy()
    {
        final HeadNode clone = new HeadNode();

        final Queue<Pair<Node, Node>> worklist = new ArrayDeque<>();
        for(Node successor : successors())
            worklist.add(Pair.of(clone, successor));

        while(!worklist.isEmpty())
        {
            final Pair<Node, Node> pair = worklist.poll();
            final Node previousNode = pair.getLeft();
            final Node toClone = pair.getRight();

            final Node newNode = new Node(toClone.Base);
            newNode.Quality = toClone.Quality;
            newNode.MaxQuality = toClone.MaxQuality;
            newNode.Support = new ArrayList<>(toClone.Support);
            previousNode.setNext(newNode);
            for(Node successor : toClone.successors())
                worklist.add(Pair.of(newNode, successor));
        }

        return clone;
    }

    private static int getReadStartIndex(final Read read, final int startPosition)
    {
        // CHECK: assumes that the read's cigar hasn't changed
        final int startIndex = ReadUtils.getReadPositionAtReferencePosition(read, startPosition) - 1;
        if(startIndex >= 0)
            return startIndex;

        int readIndex = 0;
        int refPosition = read.getAlignmentStart();

        int lastAlignedRefPosition = -1;
        int lastAlignedReadIndex = -1;

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            if(element.getOperator() == M && refPosition < startPosition)
            {
                lastAlignedRefPosition = refPosition + 1;
                lastAlignedReadIndex = readIndex + 1;
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();
        }

        if(lastAlignedRefPosition > 0)
            return (startPosition - lastAlignedRefPosition) + lastAlignedReadIndex - 1;

        return startPosition - read.getUnclippedStart();

        /*
        Alignment lastBlock = null;
        for(Alignment block : read.getAlignmentBlocks())
        {
            if(block.isMapped() && block.ReferenceStartPosition < startPosition)
                lastBlock = block;
        }

        if(lastBlock != null)
            return (startPosition - lastBlock.ReferenceStartPosition) + lastBlock.SequenceStartPosition - 1;

        return startPosition - read.getUnclippedStart();
        */
    }

    public static HeadNode createIndexed(final Read alignment, final int startIndex, final Direction orientation)
    {
        final HeadNode root = new HeadNode();
        Node current = root;
        final int startOffset = orientation == Direction.FORWARDS
                ? startIndex
                : alignment.getLength() - startIndex - 1;
        for(int i = startOffset; i < alignment.getLength(); i++)
        {
            final int index;
            if(orientation == Direction.FORWARDS)
                index = i;
            else
                index = alignment.getLength() - 1 - i;
            final byte base = alignment.getBases()[index];
            final int quality = base == 'N' ? 0 : alignment.getBaseQuality()[index];

            final Node newNode = new Node((char) (base == 'N' ? 'A' : base));
            newNode.MaxQuality = newNode.Quality = quality;
            newNode.Support = new ArrayList<>();
            newNode.Support.add(new Node.Support(alignment, index));
            current.setNext(newNode);
            current = newNode;
        }

        return root;
    }

    @Nullable
    public static HeadNode create(final Read read, final int startPosition, final boolean isForwards)
    {
        return create(read, startPosition, isForwards ? Direction.FORWARDS : Direction.REVERSE);
    }

    @Nullable
    public static HeadNode create(final Read read, final int startPosition, final Direction orientation)
    {
        final int startIndex = getReadStartIndex(read, startPosition);
        if(startIndex < 0 || startIndex >= read.getLength())
            return null;
        return createIndexed(read, startIndex, orientation);
    }

    public static HeadNode create(final SupportedAssembly assembly, final Direction orientation)
    {
        // For each base in the assembly, determine what supports it at that base
        final List<List<Node.Support>> nodeSupport = new ArrayList<>(assembly.Assembly.length());
        for(int i = 0; i < assembly.Assembly.length(); i++)
            nodeSupport.add(new ArrayList<>());

        for(ReadSupport readSupport : assembly.readSupport())
        {
            Read supportRead = readSupport.Read;
            int supportIndex = readSupport.Index;

            int length = Math.min(assembly.Assembly.length() - supportIndex, supportRead.getLength());
            if(supportIndex < 0)
            {
                length += supportIndex;
                supportIndex = 0;
            }

            for(int i = 0; i < length; i++)
                nodeSupport.get(supportIndex + i).add(new Node.Support(supportRead, i));
        }

        for(List<Node.Support> list : nodeSupport)
            if(list.size() > 1)
                list.sort(Comparator.naturalOrder());

        final HeadNode root = new HeadNode();
        Node current = root;
        for(int i = 0; i < assembly.Assembly.length(); i++)
        {
            final int index = orientation == Direction.FORWARDS ? i : assembly.Assembly.length() - 1 - i;
            final char base = assembly.Assembly.charAt(index);

            final Node newNode = new Node(base);
            newNode.Support = nodeSupport.get(i);
            newNode.recomputeQuality();
            current.setNext(newNode);
            current = newNode;
        }

        return root;
    }

    private static class MergePair
    {
        @Nullable
        public final Node Left, Right;
        public final Node Previous;

        public MergePair(@Nullable final Node left, @Nullable final Node right, final Node previous)
        {
            Left = left;
            Right = right;
            Previous = previous;
        }
    }

    /** May mutate left. Directly overlays the two nodes from position 0. */
    public static HeadNode combine(final HeadNode left, final HeadNode right)
    {
        return combine(left, right, true);
    }

    /** May mutate left. Directly overlays the two nodes from position 0. */
    public static HeadNode combine(@Nullable final HeadNode left, @Nullable final HeadNode right, final boolean sortSupport)
    {
        if(left == right || right == null)
            return left;
        else if(left == null)
            return right;

        final HeadNode result = new HeadNode();

        final Queue<MergePair> mergePairs = new ArrayDeque<>();
        mergePairs.add(new MergePair(left.nextA, right.nextA, result));
        mergePairs.add(new MergePair(left.nextT, right.nextT, result));
        mergePairs.add(new MergePair(left.nextC, right.nextC, result));
        mergePairs.add(new MergePair(left.nextG, right.nextG, result));

        while(!mergePairs.isEmpty())
        {
            final MergePair pair = mergePairs.poll();
            final Node node;
            if(pair.Left != null && pair.Right != null)
            {
                node = new Node(pair.Left.Base);
                if(pair.Left != pair.Right)
                {
                    if(sortSupport)
                    {
                        node.Support = mergeSupport(pair.Left.Support, pair.Right.Support);
                        node.recomputeQuality();
                    }
                    else
                    {
                        node.Support = pair.Left.Support instanceof ArrayList ? pair.Left.Support : new ArrayList<>(pair.Left.Support);
                        if(node.Support.size() < 10_000)
                        {
                            node.Support.addAll(pair.Right.Support);
                            if(node.Support.size() > 50_000)
                            {
                                final Set<Node.Support> supportSet = Collections.newSetFromMap(new IdentityHashMap<>());
                                supportSet.addAll(node.Support);
                                node.Support = new ArrayList<>(supportSet);
                            }
                        }
                    }

                    mergePairs.add(new MergePair(pair.Left.nextA, pair.Right.nextA, node));
                    mergePairs.add(new MergePair(pair.Left.nextT, pair.Right.nextT, node));
                    mergePairs.add(new MergePair(pair.Left.nextC, pair.Right.nextC, node));
                    mergePairs.add(new MergePair(pair.Left.nextG, pair.Right.nextG, node));
                }
                else
                {
                    node.Support = pair.Left.Support;
                    node.Quality = pair.Left.Quality;
                    node.MaxQuality = pair.Left.MaxQuality;
                    pair.Left.successors().forEach(node::setNext);
                }
            }
            else if(pair.Left != null)
            {
                left.mJunctions++;
                node = pair.Left;
            }
            else if(pair.Right != null)
            {
                right.mJunctions++;
                node = pair.Right;
            }
            else
                continue; // Both sides null
            pair.Previous.setNext(node);
        }

        return result;
    }

    public void sortSupport()
    {
        final Queue<Node> toProcess = new ArrayDeque<>();
        final Set<Node> visited = Collections.newSetFromMap(new IdentityHashMap<>());
        toProcess.add(this);

        while(!toProcess.isEmpty())
        {
            final Node node = toProcess.poll();
            if(!visited.add(node))
                continue;

            toProcess.addAll(node.successors());
            if(node instanceof HeadNode)
                continue;

            node.Support = new ArrayList<>(new LinkedHashSet<>(node.Support));
            node.Support.sort(Comparator.naturalOrder());
            node.recomputeQuality();
        }
    }
}
