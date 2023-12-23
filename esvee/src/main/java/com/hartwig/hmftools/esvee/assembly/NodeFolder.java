package com.hartwig.hmftools.esvee.assembly;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.function.BiConsumer;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.esvee.SVAConfig;
import com.hartwig.hmftools.esvee.util.Timeout;

import org.jetbrains.annotations.Nullable;

@NotThreadSafe
public class NodeFolder
{
    private final SVAConfig mConfig;
    private final Timeout mTimeout;
    private int mOperations;

    public NodeFolder(final SVAConfig config)
    {
        this(config, new Timeout(false, 0));
    }

    public NodeFolder(final SVAConfig config, final Timeout timeout)
    {
        mConfig = config;
        mTimeout = timeout;
    }

    public void prepruneNodes(final HeadNode head)
    {
        final Queue<Node> toProcess = new ArrayDeque<>();
        toProcess.add(head);

        final Set<Node> processed = Collections.newSetFromMap(new IdentityHashMap<>());
        while(!toProcess.isEmpty())
        {
            final Node node = toProcess.poll();
            if (!processed.add(node))
                continue;

            final var successors = node.successors();
            if(successors.size() == 1)
                toProcess.add(successors.get(0));
            else if(successors.size() > 1)
            {
                tryPrepruneNode(node);
                toProcess.addAll(node.successors());
            }
        }
    }

    private void tryPrepruneNode(final Node node)
    {
        final var successors = node.successors();
        if (successors.size() != 2)
            return;

        Node left = successors.get(0);
        Node right = successors.get(1);
        if (!left.successors().equals(right.successors()))
            return;

        if (left.MaxQuality > right.MaxQuality) // Ensure left has the lowest MaxQuality
        {
            final Node temp = left;
            left = right;
            right = temp;
        }

        if (left.MaxQuality > mConfig.lowBaseQualThreshold())
            return;

        if (left.Quality < right.Quality && (left.Support.size() == 1 || left.Quality < mConfig.lowBaseQualCumulativeThreshold()))
            node.removeSuccessor(left.Base);
    }

    public boolean foldPaths(final HeadNode head)
    {
        boolean madeChanges = false;
        final Queue<Node> toProcess = new ArrayDeque<>();
        toProcess.add(head);

        final Set<Node> processed = Collections.newSetFromMap(new IdentityHashMap<>());
        final Map<Node, Node> rewriteMap = new IdentityHashMap<>();

        while(!toProcess.isEmpty())
        {
            Node node = toProcess.poll();
            node = rewriteMap.getOrDefault(node, node);
            if (!processed.add(node))
                continue;

            final var successors = node.successors();
            if(successors.size() == 1)
                toProcess.add(successors.get(0));
            else if(successors.size() > 1)
            {
                madeChanges |= tryFoldJunction(node, rewriteMap);
                toProcess.addAll(node.successors());
            }
        }

        return madeChanges;
    }

    private boolean tryFoldJunction(final Node node, final Map<Node, Node> rewriteMap)
    {
        boolean madeChanges = false;
        var successors = node.successors();
        for(int i = 0; i < successors.size(); i++)
        {
            for(int j = 0; j < successors.size(); j++)
            {
                final Node left = successors.get(i);
                final Node right = successors.get(j);

                if(i == j || left.Quality < right.Quality)
                    continue; // Instead of looking at (left, right), consider (right, left).
                //noinspection SlowListContainsAll
                if (left.successors().containsAll(right.successors()))
                    continue; // Already folded

                if (tryFoldSuccessor(node, left, right, rewriteMap))
                {
                    madeChanges = true;
                    i = 0;
                    j = 0;
                    successors = node.successors();
                }
            }
        }

        return madeChanges;
    }

    private boolean tryFoldSuccessor(final Node parent, final Node left, final Node right, final Map<Node, Node> rewriteMap)
    {
        mOperations = 0;
        if (!canOverlay(left, right, 0))
            return false;

        // We will fold right's successors into left
        final List<Node> newSuccessors = new ArrayList<>();
        final BiConsumer<Node, Node> merger = (l, r) ->
        {
            if (l == r)
            {
                if (l != null)
                    newSuccessors.add(l);
                return;
            }

            final Node merged = merge(l, r);
            rewriteMap.put(l, merged);
            rewriteMap.put(r, merged);
            if (merged != null)
                newSuccessors.add(merged);
        };

        merger.accept(left.nextA, right.nextA);
        merger.accept(left.nextT, right.nextT);
        merger.accept(left.nextC, right.nextC);
        merger.accept(left.nextG, right.nextG);
        merger.accept(left.nextX, right.nextX);

        newSuccessors.forEach(left::setNext);
        for(final Node newSuccessor : newSuccessors)
            if (right.getNext(newSuccessor.Base) != null)
                right.setNext(newSuccessor);

        for(final Node sibling : parent.successors())
        {
            if (sibling == left || sibling == right)
                continue;

            final List<Node> oldCousins = sibling.successors();
            for(final Node oldCousin : oldCousins)
                if(rewriteMap.containsKey(oldCousin))
                    sibling.setNext(rewriteMap.get(oldCousin));
        }

        return true;
    }

    private boolean canOverlay(final Node left, final Node right, final int mismatchCountSoFar)
    {
        if ((mOperations++ & 0x00FFF) == 0)
            mTimeout.checkTimeout();

        if (left == null || right == null)
            return false;
        if (left == right)
            return true;

        final boolean isLowQuality = isLowQualityBase(left, true) || isLowQualityBase(right, false);
        final boolean mismatches = left.Base != right.Base && !isLowQuality;
        final int mismatchCount = mismatchCountSoFar + (mismatches ? 1 : 0);
        if (mismatchCount > mConfig.maxMismatchesForFolding())
            return false;

        // left & right can be overlayed iff every possible path from right has a candidate path through left that does not disagree
        for(final Node rightSuccessor : right.successors())
        {
            // Try the matching base first -- most likely to work
            @Nullable
            final Node candidateLeft = left.getNext(rightSuccessor.Base);
            if (canOverlay(candidateLeft, rightSuccessor, mismatchCount))
                continue;

            if (left.successors().isEmpty())
                continue; // No bases means no disagreement

            boolean couldOverlay = false;
            for (final Node leftSuccessor : left.successors())
            {
                if (leftSuccessor.Base == rightSuccessor.Base)
                    continue; // Already checked
                if (canOverlay(leftSuccessor, rightSuccessor, mismatchCount))
                {
                    couldOverlay = true;
                    break;
                }
            }
            if (!couldOverlay)
                return false;
        }

        return true;
    }

    private boolean isLowQualityBase(final Node node, final boolean checkDepth)
    {
        return node.MaxQuality <= mConfig.lowBaseQualThreshold()
                || (checkDepth && node.supportDepth() > 1 && node.Quality <= mConfig.lowBaseQualCumulativeThreshold());
    }

    @SuppressWarnings("ConstantValue")
    @Nullable
    private Node merge(@Nullable final Node left, @Nullable final Node right)
    {
        if (left == right)
            return left;
        else if(left == null && right != null)
            return right;
        else if(left != null && right == null)
            return left;
        else if(left != null && right != null)
            return fold(left, right);
        else
            return null;
    }

    Node fold(final Node destination, final Node eliminated)
    {
        final HeadNode leftHead = new HeadNode(mConfig);
        destination.successors().forEach(leftHead::setNext);

        final HeadNode rightHead = new HeadNode(mConfig);
        eliminated.successors().forEach(rightHead::setNext);

        final HeadNode merged = HeadNode.combine(leftHead, rightHead);

        final Node replacement = new Node(destination.Base);
        if (destination.Base == eliminated.Base)
        {
            replacement.Support = HeadNode.mergeSupport(destination.Support, eliminated.Support);
            replacement.recomputeQuality();
        }
        else
        {
            replacement.MaxQuality = destination.MaxQuality;
            replacement.Quality = destination.Quality;
            replacement.Support = destination.Support;
        }

        merged.successors().forEach(replacement::setNext);

        return replacement;
    }
}
