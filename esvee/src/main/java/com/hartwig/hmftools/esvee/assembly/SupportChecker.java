package com.hartwig.hmftools.esvee.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.sequence.Sequence;
import com.hartwig.hmftools.esvee.processor.SequenceDecomposer;

import org.jetbrains.annotations.Nullable;

/** Checks two sequences against each-other to see whether the right sequence supports the left sequence.
 * This check is not symmetric, as we take into account low-quality bases on the right as "neutral" bases
 * (as opposed to them indicating a disagreement) while we ignore quality scores on the left. The class
 * primarily works by decomposing sequences into runs of bases and sections that repeat (see SequenceDecomposer). */
public class SupportChecker
{
    private final int mLowBaseQualThreshold;

    public final SupportCheckerConfig StrongSupport;
    public final SupportCheckerConfig WeakSupport;
    public final SupportCheckerConfig AssemblySupport;

    public SupportChecker()
    {
        mLowBaseQualThreshold = SvConstants.LOW_BASE_QUAL_THRESHOLD;

        StrongSupport = new SupportCheckerConfig(SvConstants.MAXMISMATCHEDCOUNTFORSTRONGSUPPORT, false);
        WeakSupport = new SupportCheckerConfig(SvConstants.MAXMISMATCHEDCOUNTFORWEAKSUPPORT, false);
        AssemblySupport = new SupportCheckerConfig(SvConstants.MAXMISMATCHEDCOUNTFORDEDUPINGASSEMBLIES, true);
    }

    public class SupportCheckerConfig
    {
        private final int mEditDistance;
        private final boolean mIgnoreQual; // FIXME: Wire this through

        public SupportCheckerConfig(final int editDistance, final boolean ignoreQual)
        {
            mEditDistance = editDistance;
            mIgnoreQual = ignoreQual;
        }

        public boolean supports(final Sequence left, final Sequence right)
        {
            return supports(left, right, Math.min(left.getLength(), right.getLength()) - mEditDistance);
        }

        public boolean supports(final Sequence left, final Sequence right, final int minOverlap)
        {
            return supportIndex(left, right, minOverlap) != null;
        }

        public boolean supportsAt(final Sequence left, final Sequence right, final int index)
        {
            return supportsAtIndex(left, right, Math.min(left.getLength(), right.getLength()) - mEditDistance,
                    mEditDistance, index) > 0;
        }

        @Nullable
        public Integer supportIndex(final Sequence left, final Sequence right)
        {
            return supportIndex(left, right, Math.min(left.getLength(), right.getLength()) - mEditDistance);
        }

        @Nullable
        public Integer supportIndex(final Sequence left, final Sequence right, final int minOverlap)
        {
            return supportIndex(left, right, minOverlap, -right.getLength(), left.getLength());
        }

        @Nullable
        public Integer supportIndex(final Sequence left, final Sequence right, final int minOverlap,
                final int minCheckIndex, final int maxCheckIndex)
        {
            return determineSupportIndex(left, right, minOverlap, mEditDistance, minCheckIndex, maxCheckIndex);
        }

        public Integer bestSupportIndex(final Sequence left, final Sequence right, final int minOverlap)
        {
            return bestSupportIndex(left, right, minOverlap, -right.getLength(), left.getLength());
        }

        @Nullable
        public Integer bestSupportIndex(final Sequence left, final Sequence right, final int minOverlap,
                final int minCheckIndex, final int maxCheckIndex)
        {
            return determineBestSupportIndex(left, right, minOverlap, mEditDistance, minCheckIndex, maxCheckIndex);
        }
    }

    @Nullable
    public Integer determineSupportIndex(final Sequence left, final Sequence right, final int minOverlap, final int maxEditDistance)
    {
        return determineSupportIndex(left, right, minOverlap, maxEditDistance, -right.getLength(), left.getLength());
    }

    @Nullable
    private Integer determineSupportIndex(final Sequence left, final Sequence right, final int minOverlap, final int maxEditDistance,
            final int minCheckIndex, final int maxCheckIndex)
    {
        {
            final int startIndex = Math.max(0, minCheckIndex);
            final int stopIndex = Math.min(maxCheckIndex, left.getLength());
            for(int assemblyOffset = startIndex; assemblyOffset < stopIndex; assemblyOffset++)
                if(supportsAtIndex(left, right, minOverlap, maxEditDistance, assemblyOffset) > 0)
                    return assemblyOffset;
        }

        if(minCheckIndex < 0)
        {
            final int startIndex = Math.min(-1, maxCheckIndex);
            final int stopIndex = Math.max(minCheckIndex, -right.getLength());
            for(int assemblyOffset = startIndex; assemblyOffset > stopIndex; assemblyOffset--)
                if(supportsAtIndex(left, right, minOverlap, maxEditDistance, assemblyOffset) > 0)
                    return assemblyOffset;
        }

        return null;
    }

    @Nullable
    public Integer determineBestSupportIndex(final Sequence left, final Sequence right, final int minOverlap, final int maxEditDistance)
    {
        return determineBestSupportIndex(left, right, minOverlap, maxEditDistance, -right.getLength(), left.getLength());
    }

    @Nullable
    public Integer determineBestSupportIndex(final Sequence left, final Sequence right, final int minOverlap, final int maxEditDistance,
            final int minCheckIndex, final int maxCheckIndex)
    {
        int bestScore = 0;
        int bestAssemblyOffset = -1;
        {
            final int startIndex = Math.max(0, minCheckIndex);
            final int stopIndex = Math.min(maxCheckIndex, left.getLength());
            for(int assemblyOffset = startIndex; assemblyOffset < stopIndex; assemblyOffset++)
            {
                final int score = supportsAtIndex(left, right, minOverlap, maxEditDistance, assemblyOffset);
                if(score > bestScore)
                {
                    bestScore = score;
                    bestAssemblyOffset = assemblyOffset;
                }
            }
        }

        if(minCheckIndex < 0)
        {
            final int startIndex = Math.min(-1, maxCheckIndex);
            final int stopIndex = Math.max(minCheckIndex, -right.getLength());
            for(int assemblyOffset = startIndex; assemblyOffset > stopIndex; assemblyOffset--)
            {
                final int score = supportsAtIndex(left, right, minOverlap, maxEditDistance, assemblyOffset);
                if(score > bestScore)
                {
                    bestScore = score;
                    bestAssemblyOffset = assemblyOffset;
                }
            }
        }

        return bestScore == 0 ? null : bestAssemblyOffset;
    }

    private class SupportState
    {
        public final List<SequenceDecomposer.Node> LeftNodes;
        public final List<SequenceDecomposer.Node> RightNodes;

        public int LeftIndex, RightIndex;
        public int LeftOffset, RightOffset;

        public int AgreeBases;
        public int DisagreeBases;
        public int NeutralBases;

        public SupportState(final List<SequenceDecomposer.Node> leftNodes, final List<SequenceDecomposer.Node> rightNodes)
        {
            LeftNodes = leftNodes;
            RightNodes = rightNodes;
        }

        public SupportState(final List<SequenceDecomposer.Node> leftNodes, final List<SequenceDecomposer.Node> rightNodes,
                final SupportState other)
        {
            LeftNodes = leftNodes;
            RightNodes = rightNodes;

            LeftIndex = other.LeftIndex;
            RightIndex = other.RightIndex;
            LeftOffset = other.LeftOffset;
            RightOffset = other.RightOffset;
            AgreeBases = other.AgreeBases;
            DisagreeBases = other.DisagreeBases;
            NeutralBases = other.NeutralBases;
        }

        void advanceLeftIfRequired()
        {
            while(LeftIndex < LeftNodes.size() && LeftOffset >= LeftNodes.get(LeftIndex).length())
            {
                LeftOffset -= LeftNodes.get(LeftIndex).length();
                LeftIndex++;
            }
        }

        void advanceRightIfRequired()
        {
            while(RightIndex < RightNodes.size() && RightOffset >= RightNodes.get(RightIndex).length())
            {
                RightOffset -= RightNodes.get(RightIndex).length();
                RightIndex++;
            }
        }

        boolean isAtEnd()
        {
            return LeftIndex >= LeftNodes.size() || RightIndex >= RightNodes.size();
        }

        boolean supports(final int minOverlap)
        {
            return AgreeBases >= minOverlap && NeutralBases * 3 <= AgreeBases;
        }

        private boolean compareBasesNodes(final SequenceDecomposer.BasesNode leftNode, final SequenceDecomposer.BasesNode rightNode,
                final int lowBaseQualThreshold, final int maxEditDistance)
        {
            final byte[] leftBases = leftNode.Bases;
            final byte[] rightBases = rightNode.Bases;
            final byte[] rightQuality = rightNode.Quality;

            final int compareLength = Math.min(leftBases.length - LeftOffset, rightBases.length - RightOffset);
            for(int i = 0; i < compareLength; i++)
            {
                final byte leftBase = leftBases[i + LeftOffset];
                final byte rightBase = rightBases[i + RightOffset];
                final byte rightQual = rightQuality[i + RightOffset];

                if(leftBase == rightBase)
                    AgreeBases++;
                else
                {
                    if(rightQual <= lowBaseQualThreshold) // We only check the RHS for quality.
                        NeutralBases++;
                    else if(++DisagreeBases > maxEditDistance)
                        return false;
                }
            }

            LeftOffset += compareLength;
            RightOffset += compareLength;
            return true;
        }

        private boolean compareRepeatNodes(final SequenceDecomposer.RepeatNode leftNode, final SequenceDecomposer.RepeatNode rightNode,
                final int maxEditDistance, final boolean first)
        {
            if(LeftOffset % leftNode.Bases.length != 0 || RightOffset % rightNode.Bases.length != 0)
                return false;
            if(!Arrays.equals(leftNode.Bases, rightNode.Bases))
                return false;

            final int leftSize = leftNode.length() - LeftOffset;
            final int rightSize = rightNode.length() - RightOffset;
            final int leftRepCount = leftSize / leftNode.Bases.length;
            final int rightRepCount = rightSize / rightNode.Bases.length;
            if(first && leftRepCount != rightRepCount)
                return false;

            final int delta = Math.abs(leftRepCount - rightRepCount);
            if(delta == 0)
                AgreeBases += leftSize;
            else
            {
                AgreeBases += Math.min(leftSize, rightSize);
                // Is this at the end? If so, not a disagreement
                final boolean isLeftCutOff = LeftIndex == LeftNodes.size() - 1 && leftSize < rightSize;
                final boolean isRightCutOff = RightIndex == RightNodes.size() - 1 && leftSize > rightSize;
                if(!isLeftCutOff && !isRightCutOff)
                {
                    if(!tryProcessSplitRepeat(leftNode, rightNode, delta, maxEditDistance))
                        return false;
                }
            }

            LeftIndex++;
            RightIndex++;
            LeftOffset = RightOffset = 0;
            return true;
        }

        private boolean tryProcessSplitRepeat(final SequenceDecomposer.RepeatNode leftNode, final SequenceDecomposer.RepeatNode rightNode,
                final int delta, final int maxEditDistance)
        {
            if(delta > 1)
            {
                if(leftNode.RepeatCount > rightNode.RepeatCount)
                {
                    // Peek right
                    final var nextRight = peekRight(1);
                    assert nextRight != null : "Should not have been called without a least 1 successor node";

                    if(nextRight instanceof SequenceDecomposer.BasesNode)
                    {

                    }
                }
                else
                {
                    // Peek left

                }
            }

            final int differenceScore = repeatDifferenceScore(leftNode.RepeatCount, rightNode.RepeatCount);
            DisagreeBases += differenceScore;
            return DisagreeBases <= maxEditDistance;
        }

        private int repeatDifferenceScore(final int leftLength, final int rightLength)
        {
            // TODO: Consider penalising long repeat differences less than short ones (eg 5 -> 7 is worse than 14 -> 16)
            final int delta = Math.abs(leftLength - rightLength);
            return delta / 2 + 1;
        }

        @Nullable
        private SequenceDecomposer.Node peekLeft(final int peek)
        {
            final int peekIndex = LeftIndex + peek;
            if(peekIndex >= LeftNodes.size())
                return null;
            return LeftNodes.get(peekIndex);
        }

        @Nullable
        private SequenceDecomposer.Node peekRight(final int peek)
        {
            final int peekIndex = RightIndex + peek;
            if(peekIndex >= RightNodes.size())
                return null;
            return RightNodes.get(peekIndex);
        }

        private SupportState pushLeft(final List<SequenceDecomposer.Node> nodes, final boolean removeTwo)
        {
            final List<SequenceDecomposer.Node> leftNodes = new ArrayList<>(LeftNodes);
            leftNodes.addAll(LeftIndex, nodes);
            leftNodes.remove(LeftIndex + nodes.size());
            if(removeTwo)
                leftNodes.remove(LeftIndex + nodes.size());

            final SupportState newState = new SupportState(leftNodes, RightNodes, this);
            newState.LeftOffset = 0;
            return newState;
        }

        private SupportState pushRight(final List<SequenceDecomposer.Node> nodes, final boolean removeTwo)
        {
            final List<SequenceDecomposer.Node> rightNodes = new ArrayList<>(RightNodes);
            rightNodes.addAll(RightIndex, nodes);
            rightNodes.remove(RightIndex + nodes.size());
            if(removeTwo)
                rightNodes.remove(RightIndex + nodes.size());

            final SupportState newState = new SupportState(LeftNodes, rightNodes, this);
            newState.RightOffset = 0;
            return newState;
        }

        @Nullable
        public SequenceDecomposer.BasesNode nextLeftBases()
        {
            if(LeftIndex + 1 >= LeftNodes.size())
                return null;
            final SequenceDecomposer.Node next = LeftNodes.get(LeftIndex + 1);
            return next instanceof SequenceDecomposer.BasesNode ? (SequenceDecomposer.BasesNode) next : null;
        }

        @Nullable
        public SequenceDecomposer.BasesNode nextRightBases()
        {
            if(RightIndex + 1 >= RightNodes.size())
                return null;
            final SequenceDecomposer.Node next = RightNodes.get(RightIndex + 1);
            return next instanceof SequenceDecomposer.BasesNode ? (SequenceDecomposer.BasesNode) next : null;
        }
    }

    public int supportsAtIndex(final Sequence left, final Sequence right, final int minOverlap, final int maxEditDistance,
            final int checkIndex)
    {
        SupportState state = new SupportState(left.decompose(), right.decompose());
        state.LeftOffset = checkIndex;
        if(checkIndex < 0)
        {
            state.RightOffset = -checkIndex;
            state.LeftOffset = 0;
        }

        boolean first = true;
        while(true)
        {
            state.advanceLeftIfRequired();
            state.advanceRightIfRequired();
            if(state.isAtEnd())
                break;

            final SequenceDecomposer.Node leftNode = state.LeftNodes.get(state.LeftIndex);
            final SequenceDecomposer.Node rightNode = state.RightNodes.get(state.RightIndex);
            if(leftNode instanceof SequenceDecomposer.BasesNode && rightNode instanceof SequenceDecomposer.BasesNode)
            {
                if(!state.compareBasesNodes((SequenceDecomposer.BasesNode) leftNode, (SequenceDecomposer.BasesNode) rightNode,
                        mLowBaseQualThreshold, maxEditDistance))
                    return 0;
            }
            else if(leftNode instanceof SequenceDecomposer.RepeatNode && rightNode instanceof SequenceDecomposer.RepeatNode)
            {
                final SequenceDecomposer.RepeatNode leftRepeat = (SequenceDecomposer.RepeatNode) leftNode;
                final SequenceDecomposer.RepeatNode rightRepeat = (SequenceDecomposer.RepeatNode) rightNode;
                if(state.LeftOffset % leftRepeat.Bases.length != 0) // Left is out-of-phase
                {
                    @Nullable
                    final SequenceDecomposer.BasesNode next = state.nextLeftBases();
                    final List<SequenceDecomposer.Node> shifted = leftRepeat.shift(state.LeftOffset, 0, next);
                    state = state.pushLeft(shifted, next != null);
                    continue;
                }
                if(state.RightOffset % rightRepeat.Bases.length != 0) // Right is out-of-phase
                {
                    @Nullable
                    final SequenceDecomposer.BasesNode next = state.nextRightBases();
                    final List<SequenceDecomposer.Node> shifted = rightRepeat.shift(state.RightOffset, 0, next);
                    state = state.pushRight(shifted, next != null);
                    continue;
                }
                if(!state.compareRepeatNodes(leftRepeat, rightRepeat, maxEditDistance, first))
                    return 0;
            }
            else
            {
                final boolean leftIsBases = leftNode instanceof SequenceDecomposer.BasesNode;
                // One side is repeat, the other is bases
                //noinspection DataFlowIssue
                final SequenceDecomposer.RepeatNode repeatNode = leftIsBases
                        ? (SequenceDecomposer.RepeatNode) rightNode
                        : (SequenceDecomposer.RepeatNode) leftNode;
                final SequenceDecomposer.BasesNode basesNode = leftIsBases
                        ? (SequenceDecomposer.BasesNode) leftNode
                        : (SequenceDecomposer.BasesNode) rightNode;

                final int basesOffset = leftIsBases ? state.LeftOffset : state.RightOffset;
                final int repeatOffset = leftIsBases ? state.RightOffset : state.LeftOffset;
                // Decompose the repeat node by up to the bases length
                //noinspection DataFlowIssue
                final int size = Math.min(repeatNode.length() - repeatOffset, basesNode.length() - basesOffset);
                @Nullable
                final SequenceDecomposer.BasesNode next = leftIsBases ? state.nextRightBases() : state.nextLeftBases();
                final List<SequenceDecomposer.Node> shifted = repeatNode.shift(repeatOffset, size, next);
                if(leftIsBases)
                    state = state.pushRight(shifted, next != null);
                else
                    state = state.pushLeft(shifted, next != null);
            }
            first = false;
        }

        return state.supports(minOverlap) ? (state.AgreeBases - state.DisagreeBases * 5) : 0;
    }
}
