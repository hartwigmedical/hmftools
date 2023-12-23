package com.hartwig.hmftools.esvee.processor;

import java.util.List;

import com.hartwig.hmftools.esvee.models.Sequence;

import org.apache.commons.lang3.tuple.Pair;

public class SequenceMerger
{
    private int mOutputPosition;
    private final byte[] mOutputBases;
    private final byte[] mOutputQual;

    private SequenceMerger(final int maxOutputLength)
    {
        mOutputBases = new byte[maxOutputLength];
        mOutputQual = new byte[maxOutputLength];
    }

    private Sequence mergeInner(final Sequence left, final Sequence right, final int supportIndex)
    {
        final int centerSize = Math.min(left.getLength() - supportIndex, right.getLength());

        final var leftNodes = left.decompose();
        final var rightNodes = right.decompose();

        NodePosition leftPosition = new NodePosition(leftNodes);
        NodePosition rightPosition = new NodePosition(rightNodes);

        leftPosition = copyBases(leftPosition, supportIndex);

        final var result = mergeBases(leftPosition, rightPosition, centerSize);
        leftPosition = result.getLeft();
        rightPosition = result.getRight();

        copyBases(leftPosition, Integer.MAX_VALUE);
        copyBases(rightPosition, Integer.MAX_VALUE);

        return Sequence.fromBytes(mOutputBases, mOutputQual, mOutputPosition);
    }

    /**
     * @return a SimpleSequence containing the merged product of the two sequences.
     */
    public static Sequence merge(final Sequence left, final Sequence right, final int supportIndex)
    {
        if(supportIndex < 0)
            return merge(right, left, -supportIndex);

        return new SequenceMerger(left.getLength() + right.getLength())
                .mergeInner(left, right, supportIndex);
    }

    private Pair<NodePosition, NodePosition> mergeBases(
            NodePosition leftPosition, NodePosition rightPosition,
            final int length)
    {
        int copiedBases = 0;
        while (copiedBases < length && leftPosition.hasData() && rightPosition.hasData())
        {
            final SequenceDecomposer.Node l = leftPosition.node();
            final SequenceDecomposer.Node r = rightPosition.node();

            if (l instanceof SequenceDecomposer.RepeatNode && r instanceof SequenceDecomposer.RepeatNode)
            {
                final SequenceDecomposer.RepeatNode rl = (SequenceDecomposer.RepeatNode) l;
                final SequenceDecomposer.RepeatNode rr = (SequenceDecomposer.RepeatNode) r;
                final SequenceDecomposer.RepeatNode repeat = rl.quality() > rr.quality()
                        ? rl
                        : rr;

                final int satelliteLength = repeat.Bases.length;
                for(int i = 0; i < repeat.RepeatCount; i++)
                {
                    for(int j = 0; j < satelliteLength; j++)
                    {
                        final byte base = repeat.base(j);
                        final byte qual = (byte) Math.max(rl.qual(satelliteLength * j + i), rr.qual(satelliteLength * j + i));

                        mOutputBases[mOutputPosition + j] = base;
                        mOutputQual[mOutputPosition + j] = qual;
                    }
                    mOutputPosition += satelliteLength;
                }
                copiedBases += (rl.RepeatCount * satelliteLength);
                leftPosition = leftPosition.advanceNode();
                rightPosition = rightPosition.advanceNode();
            }
            else
            {
                final int remainingLeftBases = l.length() - leftPosition.Offset;
                final int remainingRightBases = r.length() - rightPosition.Offset;
                final int consumeBases = Math.min(Math.min(remainingLeftBases, remainingRightBases), length - copiedBases);

                for(int i = 0; i < consumeBases; i++)
                {
                    if (l.qual(i) > r.qual(i))
                    {
                        mOutputBases[mOutputPosition + i] = l.base(leftPosition.Offset + i);
                        mOutputQual[mOutputPosition + i] = l.qual(leftPosition.Offset + i);
                    }
                    else
                    {
                        mOutputBases[mOutputPosition + i] = r.base(rightPosition.Offset + i);
                        mOutputQual[mOutputPosition + i] = r.qual(rightPosition.Offset + i);
                    }
                }
                mOutputPosition += consumeBases;

                copiedBases += consumeBases;
                leftPosition = leftPosition.advanceOffset(consumeBases);
                rightPosition = rightPosition.advanceOffset(consumeBases);
            }
        }

        return Pair.of(leftPosition, rightPosition);
    }

    private NodePosition copyBases(NodePosition position, final int length)
    {
        int copiedBases = 0;
        while (copiedBases < length && position.hasData())
        {
            final SequenceDecomposer.Node node = position.node();
            final int nodeBases = node.length() - position.Offset;
            final int consumeBases = Math.min(nodeBases, length - copiedBases);

            for(int i = 0; i < consumeBases; i++)
            {
                mOutputBases[mOutputPosition + i] = node.base(position.Offset + i);
                mOutputQual[mOutputPosition + i] = node.qual(position.Offset + i);
            }
            mOutputPosition += consumeBases;

            copiedBases += consumeBases;
            position = position.advanceOffset(consumeBases);
        }
        return position;
    }

    private static class NodePosition
    {
        private final List<SequenceDecomposer.Node> Nodes;
        private final int Index;
        private final int Offset;

        private NodePosition(final List<SequenceDecomposer.Node> nodes, final int index, final int offset)
        {
            Nodes = nodes;
            Index = index;
            Offset = offset;
        }

        public NodePosition(final List<SequenceDecomposer.Node> nodes)
        {
            this(nodes, 0, 0);
        }

        public SequenceDecomposer.Node node()
        {
            return Nodes.get(Index);
        }

        public boolean hasData()
        {
            return Index < Nodes.size() && Offset < Nodes.get(Index).length();
        }

        // @Contract(pure = true)
        public NodePosition advanceOffset(final int amount)
        {
            if (Offset + amount >= Nodes.get(Index).length())
                return new NodePosition(Nodes, Index + 1, Offset + amount - Nodes.get(Index).length());
            else
                return new NodePosition(Nodes, Index, Offset + amount);
        }

        // @Contract(pure = true)
        public NodePosition advanceNode()
        {
            return new NodePosition(Nodes, Index + 1, 0);
        }
    }
}
