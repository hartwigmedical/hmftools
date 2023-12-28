package com.hartwig.hmftools.esvee.processor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.sequence.Sequence;
import com.hartwig.hmftools.esvee.sequence.SupportedAssembly;

import org.jetbrains.annotations.Nullable;

public enum SequenceDecomposer
{
    ;

    public static final int MIN_REPETITIONS_TO_CALL_REPEAT = 5;
    public static final int MAX_MICROSATELLITE_LENGTH = 6;

    public abstract static class Node
    {
        public abstract int length();
        public abstract byte base(int i);
        public abstract byte qual(int i);
    }

    public static class BasesNode extends Node
    {
        public final byte[] Bases;
        public final byte[] Quality;
        // We don't bother tracking support for base runs, as it could be misleading.

        public BasesNode(final byte[] bases, final byte[] quality)
        {
            Bases = bases;
            Quality = quality;
        }

        @Override
        public int length()
        {
            return Bases.length;
        }

        @Override
        public byte base(final int i)
        {
            return Bases[i];
        }

        @Override
        public byte qual(final int i)
        {
            return Quality[i];
        }

        @Override
        public boolean equals(final Object obj)
        {
            if(!(obj instanceof BasesNode))
                return false;

            final BasesNode other = (BasesNode) obj;
            return Arrays.equals(Bases, other.Bases);
        }

        @Override
        public int hashCode()
        {
            return Arrays.hashCode(Bases);
        }

        @Override
        public String toString()
        {
            return new String(Bases);
        }
    }

    public static class RepeatNode extends Node
    {
        public final byte[] Bases;
        public final byte[] Qual;
        public final int RepeatCount;
        public final int SupportDepth;
        private int mQualityScore = -1;

        public RepeatNode(final byte[] bases, final byte[] qual, final int repeatCount, final int supportDepth)
        {
            assert qual.length == bases.length * repeatCount;

            Bases = bases;
            Qual = qual;
            RepeatCount = repeatCount;
            SupportDepth = supportDepth;
        }

        @Override
        public byte base(final int i)
        {
            return Bases[i % Bases.length];
        }

        @Override
        public byte qual(final int i)
        {
            if(i >= Qual.length)
                return 0;
            return Qual[i];
        }

        public int quality()
        {
            if(mQualityScore != -1)
                return mQualityScore;

            int sum = 0;
            for(byte b : Qual)
                sum += b;
            final int averageQual = sum / Qual.length;

            mQualityScore = averageQual * SupportDepth;
            return mQualityScore;
        }

        @Override
        public int length()
        {
            return Bases.length * RepeatCount;
        }

        @Override
        public boolean equals(final Object obj)
        {
            if(!(obj instanceof RepeatNode))
                return false;

            final RepeatNode other = (RepeatNode) obj;

            return Arrays.equals(Bases, other.Bases) && RepeatCount == other.RepeatCount;
        }

        @Override
        public int hashCode()
        {
            return Arrays.hashCode(Bases) ^ RepeatCount;
        }

        @Override
        public String toString()
        {
            return new String(Bases) + "x" + RepeatCount;
        }

        public List<Node> shift(final int skip, final int bases, @Nullable final BasesNode next)
        {
            final byte[] leftBases = new byte[bases];
            final byte[] leftQuals = new byte[bases];
            for(int i = 0; i < bases; i++)
            {
                leftBases[i] = base(i + skip);
                leftQuals[i] = qual(i + skip);
            }
            final BasesNode left = new BasesNode(leftBases, leftQuals);
            final int usedBases = skip + bases;

            if(usedBases % Bases.length == 0)
            {
                // Exact reduction
                final int newRepeatCount = RepeatCount - (usedBases / Bases.length);
                if(newRepeatCount <= 0)
                    return nonNullNonEmpty(left, next);

                final byte[] newQual = Arrays.copyOfRange(Qual, usedBases, Qual.length);
                final RepeatNode newRepeat = new RepeatNode(Bases, newQual, newRepeatCount, SupportDepth);

                return nonNullNonEmpty(left, newRepeat, next);
            }
            else
            {
                // Left-over bases
                // We need to 're-phase' the repeat
                // ATATAT -> A, TATA, T
                final int shift = usedBases % Bases.length;
                final byte[] newBases = new byte[Bases.length];
                for(int i = 0; i < Bases.length; i++)
                    newBases[i] = base(i + shift);

                final int newRepeatCount = RepeatCount - (usedBases / Bases.length) - 1;
                if(newRepeatCount < 0)
                    return nonNullNonEmpty(left, next);

                final int leftoverCount = Bases.length - shift;
                final byte[] newQual = Arrays.copyOfRange(Qual, usedBases, Qual.length - leftoverCount);
                RepeatNode newRepeat = new RepeatNode(newBases, newQual, newRepeatCount, SupportDepth);

                BasesNode right;
                if(next == null)
                {
                    final byte[] rightBases = Arrays.copyOfRange(Bases, shift, Bases.length);
                    final byte[] rightQual = Arrays.copyOfRange(Qual, Qual.length - leftoverCount, Qual.length);
                    right = new BasesNode(rightBases, rightQual);
                }
                else
                {
                    final byte[] rightBases = new byte[leftoverCount + next.length()];
                    final byte[] rightQual = new byte[leftoverCount + next.length()];

                    System.arraycopy(Bases, shift, rightBases, 0, Bases.length - shift);
                    System.arraycopy(Qual, Qual.length - leftoverCount, rightQual, 0, Bases.length - shift);

                    System.arraycopy(next.Bases, 0, rightBases, leftoverCount, next.length());
                    System.arraycopy(next.Quality, 0, rightQual, leftoverCount, next.length());

                    right = new BasesNode(rightBases, rightQual);
                }

                while(true)
                {
                    final int repeatSequenceSize = newRepeat.Bases.length;
                    if(right.Bases.length < repeatSequenceSize ||
                            !Arrays.equals(right.Bases, 0, repeatSequenceSize, newRepeat.Bases, 0, repeatSequenceSize))
                        break;

                    // Pull a repeat out of right
                    final byte[] qual = Arrays.copyOf(newRepeat.Qual, newRepeat.Qual.length + repeatSequenceSize);
                    System.arraycopy(right.Quality, 0, qual, newRepeat.Qual.length, repeatSequenceSize);

                    final byte[] rightBases = Arrays.copyOfRange(right.Bases, repeatSequenceSize, right.Bases.length);
                    final byte[] rightQual = Arrays.copyOfRange(right.Quality, repeatSequenceSize, right.Quality.length);
                    right = new BasesNode(rightBases, rightQual);
                    newRepeat = new RepeatNode(newRepeat.Bases, qual, newRepeat.RepeatCount + 1, newRepeat.SupportDepth);
                }

                if(newRepeatCount == 0)
                    return nonNullNonEmpty(left, right);
                else
                    return nonNullNonEmpty(left, newRepeat, right);
            }
        }

        private List<Node> nonNullNonEmpty(final Node... nodes)
        {
            final List<Node> replacements = new ArrayList<>();
            for(@Nullable final Node node : nodes)
                if(node != null && node.length() > 0)
                    replacements.add(node);
            return replacements;
        }
    }

    static class BasesNodeBuilder
    {
        private byte[] mBases = new byte[100];
        private byte[] mQuals = new byte[100];
        private int mLength = 0;

        public void append(final byte base, final byte qual)
        {
            if(mLength + 1 >= mBases.length)
            {
                mBases = Arrays.copyOf(mBases, mBases.length * 2);
                mQuals = Arrays.copyOf(mQuals, mQuals.length * 2);
            }

            mBases[mLength] = base;
            mQuals[mLength] = qual;
            mLength++;
        }

        public void append(final byte[] bases, final byte[] quals, final int fromIndex, final int length)
        {
            for(int i = 0; i < length; i++)
                append(bases[fromIndex + i], quals[fromIndex + i]);
        }

        public int length()
        {
            return mLength;
        }

        public void undo(final int count)
        {
            mLength -= count;
        }

        public byte previous(final int lookback)
        {
            return lookback >= mLength ? 0 : mBases[mLength - lookback - 1];
        }

        public byte[] extract(final int lookback, final int length)
        {
            final byte[] bytes = new byte[length];
            System.arraycopy(mBases, mLength - lookback - 1, bytes, 0, length);
            return bytes;
        }

        public BasesNode buildAndReset()
        {
            final BasesNode node = new BasesNode(
                    Arrays.copyOf(mBases, mLength),
                    Arrays.copyOf(mQuals, mLength)
            );
            mLength = 0;
            return node;
        }
    }

    public static List<Node> decompose(final SupportedAssembly assembly)
    {
        return annotateDepth(decompose((Sequence) assembly), assembly);
    }

    private static List<Node> annotateDepth(final List<Node> nodes, final SupportedAssembly assembly)
    {
        int endIndex = 0;
        for(int i = 0; i < nodes.size(); i++)
        {
            final Node node = nodes.get(i);
            endIndex += node.length();
            if(!(node instanceof RepeatNode))
                continue;

            final RepeatNode repeat = (RepeatNode) node;
            final int startIndex = endIndex - node.length();

            int supportDepth = 0;
            for(Map.Entry<Read, Integer> entry : assembly.getSupport())
            {
                final int supportStartIndex = entry.getValue();
                final int supportEndIndex = supportStartIndex + entry.getKey().getLength();

                if(supportStartIndex <= startIndex && supportEndIndex >= endIndex)
                    supportDepth++;
            }

            nodes.set(i, new RepeatNode(repeat.Bases, repeat.Qual, repeat.RepeatCount, supportDepth));
        }

        return nodes;
    }

    public static List<Node> decompose(final Sequence sequence)
    {
        return decompose(sequence.getBases(), sequence.getBaseQuality());
    }

    public static List<Node> decompose(final byte[] bases, final byte[] baseQuality)
    {
        return fixLongSatellites(decomposeInner(bases, baseQuality), baseQuality);
    }

    private static List<Node> decomposeInner(final byte[] bases, final byte[] baseQuality)
    {
        final List<Node> nodes = new ArrayList<>();

        boolean inMultiBaseRepeat = false;
        int multiBaseRepeatLength = 0;

        final BasesNodeBuilder builder = new BasesNodeBuilder();

        builder.append(bases[0], baseQuality[0]);
        int runLength = 1;
        for(int i = 1; i < bases.length; i++)
        {
            final byte base = bases[i];
            final byte qual = baseQuality[i];
            if(inMultiBaseRepeat)
            {
                if(base == builder.previous(multiBaseRepeatLength - 1))
                {
                    runLength++;
                    builder.append(base, qual);
                    continue;
                }

                // C-c-c-combo breaker
                final int repeats = runLength / multiBaseRepeatLength; // Round down, TATA = 2, TATAT = 2 (with a left-over T)
                final byte[] repeatBases = builder.extract(runLength - 1, multiBaseRepeatLength);

                builder.undo(runLength);
                if(builder.length() > 0)
                    nodes.add(builder.buildAndReset());

                final int putBackBases = runLength % multiBaseRepeatLength;
                final int startIndex = i - putBackBases - (repeatBases.length * repeats);
                final byte[] repeatQual = Arrays.copyOfRange(baseQuality, startIndex, startIndex + (repeats * multiBaseRepeatLength));

                nodes.add(new RepeatNode(repeatBases, repeatQual, repeats, 1));

                // Put the current (mismatched) base, plus any partially completed copies of the run

                builder.append(bases, baseQuality, i - putBackBases, putBackBases);
                builder.append(base, qual);

                runLength = 1;
                // We may be in a single base repeat right now -- that's the only length possible given maxLen=6, minRepeats=4.
                // Increase runLength as appropriate
                for(int j = 1; j < MAX_MICROSATELLITE_LENGTH; j++)
                    if(builder.previous(j) == builder.previous(j - 1))
                        runLength++;
                    else
                        break;
                inMultiBaseRepeat = runLength >= MIN_REPETITIONS_TO_CALL_REPEAT;
                if(inMultiBaseRepeat)
                    multiBaseRepeatLength = 1;
                continue;
            }

            builder.append(base, qual);

            // Should we start a multi-base repeat?
            for(int repeatLength = 1; repeatLength <= MAX_MICROSATELLITE_LENGTH; repeatLength++)
            {
                if(builder.length() < repeatLength * MIN_REPETITIONS_TO_CALL_REPEAT)
                    break;

                if(tryStartMultiBaseRepeat(builder, repeatLength, MIN_REPETITIONS_TO_CALL_REPEAT))
                {
                    inMultiBaseRepeat = true;
                    runLength = repeatLength * MIN_REPETITIONS_TO_CALL_REPEAT;
                    multiBaseRepeatLength = repeatLength;
                }
            }
        }

        if(builder.length() > 0)
        {
            if(inMultiBaseRepeat)
            {
                // What is our repeat?
                final int repeats = runLength / multiBaseRepeatLength; // Round down, TATA = 2, TATAT = 2 (with a left-over T)
                final byte[] repeatBases = builder.extract(runLength - 1, multiBaseRepeatLength);

                // Wind back the repeat, and add in any bases that came before the repeat.
                builder.undo(runLength);
                if(builder.length() > 0)
                    nodes.add(builder.buildAndReset());

                final int putBackBases = runLength % multiBaseRepeatLength;
                final int startIndex = bases.length - putBackBases - (repeatBases.length * repeats);
                final byte[] repeatQual = Arrays.copyOfRange(baseQuality, startIndex, startIndex + (repeats * multiBaseRepeatLength));
                nodes.add(new RepeatNode(repeatBases, repeatQual, repeats, 1));

                // Copy the incomplete repetitions
                builder.append(bases, baseQuality, bases.length - putBackBases, putBackBases);
                if(builder.length() > 0)
                    nodes.add(builder.buildAndReset());
            }
            else
                nodes.add(builder.buildAndReset());
        }

        return nodes;
    }

    private static List<Node> fixLongSatellites(final List<Node> nodes, final byte[] baseQual)
    {
        //noinspection ConstantValue
        if(MIN_REPETITIONS_TO_CALL_REPEAT >= MAX_MICROSATELLITE_LENGTH)
            return nodes;

        // If MIN_REPETITIONS_TO_CALL_REPEAT < MAX_MICROSATELLITE_LENGTH, it's possible for us to mis-call repeats like
        // AAAAATAAAAATAAAAATAAAAAT as Ax5, T, Ax5, T, ..., when it should instead be AAAAATx4.
        // Rather than making the decomposition process more complex, we just look to identify when this has occurred and "correct" it.

        if(nodes.size() < MIN_REPETITIONS_TO_CALL_REPEAT * 2)
            return nodes;

        // We look for one less, then attempt to slide the repeat around in case the edge is "hiding" in a neighbouring BasesNode
        final int desiredReps = MIN_REPETITIONS_TO_CALL_REPEAT;
        int repeatStartIndex = -1;
        for(int i = desiredReps * 2 - 1; i < nodes.size(); i++)
        {
            boolean matches = true;
            for(int reps = 1; reps < desiredReps; reps++)
            {
                if(!nodes.get(i).equals(nodes.get(i - (reps * 2)))
                        || !nodes.get(i - 1).equals(nodes.get(i - 1 - (reps * 2))))
                {
                    matches = false;
                    break;
                }
            }
            if(matches)
            {
                repeatStartIndex = i - desiredReps * 2 + 1;
                break;
            }
        }
        if(repeatStartIndex == -1)
            return nodes;

        final List<Node> newNodes = new ArrayList<>();
        int currentIndex = 0;
        for(int i = 0; i < repeatStartIndex; i++)
        {
            newNodes.add(nodes.get(i));
            currentIndex += nodes.get(i).length();
        }

        final String repeatSequence = getBases(nodes.get(repeatStartIndex)) + getBases(nodes.get(repeatStartIndex + 1));
        int repeatCount = 1;
        for(int i = repeatStartIndex + 3; i < nodes.size(); i += 2)
        {
            final Node first = nodes.get(i - 1);
            final Node second = nodes.get(i);
            if(!first.equals(nodes.get(i - 3)) || !second.equals(nodes.get(i - 2)))
                break;
            repeatCount++;
        }

        final byte[] repeatQual = Arrays.copyOfRange(baseQual, currentIndex, currentIndex + (repeatCount * repeatSequence.length()));
        newNodes.add(new RepeatNode(repeatSequence.getBytes(), repeatQual, repeatCount, 1));
        for(int i = repeatStartIndex + (repeatCount * 2); i < nodes.size(); i++)
            newNodes.add(nodes.get(i));

        return fixLongSatellites(newNodes, baseQual);
    }

    private static String getBases(final Node node)
    {
        if(node instanceof BasesNode)
            return node.toString();
        final RepeatNode repeatNode = (RepeatNode) node;
        return new String(repeatNode.Bases).repeat(repeatNode.RepeatCount);
    }

    private static boolean tryStartMultiBaseRepeat(final BasesNodeBuilder builder, final int length, final int minRepetitions)
    {
        for(int j = 0; j < minRepetitions; j++)
            for(int i = 0; i < length; i++)
                if(builder.previous(i) != builder.previous(j * length + i))
                    return false;
        return true;
    }
}
