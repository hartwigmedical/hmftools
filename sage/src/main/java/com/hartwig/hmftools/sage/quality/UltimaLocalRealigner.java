package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.apache.commons.lang3.tuple.Pair;
import org.pcollections.PVector;
import org.pcollections.TreePVector;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// TODO: Integrate into classes where appropriate
// TODO: Guard calling this against config
// TODO: clean up unneeded functions.
public class UltimaLocalRealigner
{
    @VisibleForTesting
    public static class Homopolymer
    {
        // TODO: Just use byte instead of char.
        public final char Base;
        public final int Length;

        public Homopolymer(final char base, final int length)
        {
            Base = base;
            Length = length;
        }

        @Override
        public String toString()
        {
            return String.valueOf(Base) + "x" + Length;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof Homopolymer))
            {
                return false;
            }
            final Homopolymer that = (Homopolymer) o;
            return Base == that.Base && Length == that.Length;
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(Base, Length);
        }
    }

    @VisibleForTesting
    public static class HomopolymerPair
    {
        @Nullable
        public final Homopolymer RefHomopolymer;
        @Nullable
        public final Homopolymer ReadHomopolymer;

        public HomopolymerPair(@Nullable final Homopolymer refHomopolymer, @Nullable final Homopolymer readHomopolymer)
        {
            // TODO:
            assert refHomopolymer != null || readHomopolymer != null;
            assert refHomopolymer == null || readHomopolymer == null || refHomopolymer.Base == readHomopolymer.Base;

            RefHomopolymer = refHomopolymer;
            ReadHomopolymer = readHomopolymer;
        }

        public char base()
        {
            return RefHomopolymer == null ? ReadHomopolymer.Base : RefHomopolymer.Base;
        }

        public int indelLength()
        {
            if(RefHomopolymer != null && ReadHomopolymer != null)
            {
                return ReadHomopolymer.Length - RefHomopolymer.Length;
            }

            return RefHomopolymer == null ? ReadHomopolymer.Length : -RefHomopolymer.Length;
        }

        @Override
        public String toString()
        {
            String refString = RefHomopolymer == null ? "null" : RefHomopolymer.toString();
            String readString = ReadHomopolymer == null ? "null" : ReadHomopolymer.toString();
            return format("ref(%s) read(%s)", refString, readString);
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof HomopolymerPair))
            {
                return false;
            }
            final HomopolymerPair that = (HomopolymerPair) o;
            return Objects.equals(RefHomopolymer, that.RefHomopolymer) && Objects.equals(ReadHomopolymer, that.ReadHomopolymer);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(RefHomopolymer, ReadHomopolymer);
        }
    }

    // TODO: Avoid strings?
    private static List<Homopolymer> getHomopolymers(final String bases)
    {
        List<Homopolymer> homopolymers = Lists.newArrayList();
        if(bases.length() == 0)
        {
            return homopolymers;
        }

        char currentBase = bases.charAt(0);
        int currentLength = 1;
        for(int i = 1; i < bases.length(); i++)
        {
            char base = bases.charAt(i);
            if(base == currentBase)
            {
                currentLength++;
                continue;
            }

            homopolymers.add(new Homopolymer(currentBase, currentLength));
            currentBase = base;
            currentLength = 1;
        }

        homopolymers.add(new Homopolymer(currentBase, currentLength));
        return homopolymers;
    }

    // TODO: Does this ever return false? Test this.
    public static boolean isVariantExplained(final VariantReadContext readContext)
    {
        // TODO: consider indel variants.
        final SimpleVariant variant = readContext.variant();

        // TODO: Test longer dels
        // TODO: Test inserts
        assert variant.indelLength() <= 0;

        // TODO: consider mnv variants.
        assert !variant.isMNV();

        // TODO: Consider sandwiched SNV/MNV.
        assert !isSandwichedNonIndel(readContext, getVariantRefIndex(readContext));

        // TODO: Extend so that read base count matches core length modulo indel size, plus one base padding.
        // TODO: Extend read bases and core to cover full homopolymers.
        // TODO: Mask sandwich snvs/mvns.

        List<Homopolymer> refCoreHomopolymers = getHomopolymers(readContext.refBases());
        List<Homopolymer> readCoreHomopolymers = getHomopolymers(readContext.coreStr());

        // Construct read cigar in core
        String readCigar = readContext.readCigar();
        List<CigarOperator> cigarElements = cigarStringToOps(readCigar);
        List<CigarOperator> coreCigarOps = getCoreCigarOps(readContext, cigarElements);
//        List<CigarElement> coreCigarElements = collapseCigarOps(coreCigarOps);

        // TODO: Check core length.
        int actualCoreLength = readContext.coreLength();
        int expectedCoreLength = coreCigarOps.stream().mapToInt(op -> op.consumesReadBases() ? 1 : 0).sum();
        assert actualCoreLength == expectedCoreLength;

        // TODO: Check ref bases length.
        int actualRefBasesLength = readContext.RefBases.length;
        int expectedRefBasesLength = coreCigarOps.stream().mapToInt(op -> op.consumesReferenceBases() ? 1 : 0).sum();
        assert actualRefBasesLength == expectedRefBasesLength;

        // pair homopolymers between ref and read
        List<List<HomopolymerPair>> homopolymerPairs = pairHomopolymers(refCoreHomopolymers, readCoreHomopolymers);

        // TODO: Understand ultima qual calculator. Compute qual. Only if variant still exists after re-alignment.
        return homopolymerPairs.stream()
                .map(pairs -> isVariantExplainedHelper(readContext, pairs))
                .reduce(false, (acc, x) -> acc || x);
    }

    // TODO: come up with a better way
    // TODO: Really need to come up with a better method.
    private static boolean isVariantExplainedHelper(final VariantReadContext readContext, final List<HomopolymerPair> homopolymerPairs)
    {
        SimpleVariant variant = readContext.variant();

        // TODO:
        assert variant.indelLength() <= 0;
        assert !variant.isMNV();

        // now check if homopolymers indels explain the variant.
        if(variant.isSNV())
        {
            // get index of homopolymer pair associated to variant.
            int varCoreIndex = readContext.VarIndex - readContext.CoreIndexStart + 1;
            int varHomopolymerPairIndex = -1;

            int readBasesConsumed = 0;
            while(readBasesConsumed < varCoreIndex)
            {
                Homopolymer readHomopolymer = homopolymerPairs.get(++varHomopolymerPairIndex).ReadHomopolymer;
                if(readHomopolymer != null)
                {
                    readBasesConsumed += readHomopolymer.Length;
                }
            }

            HomopolymerPair varHomopolymerPair = homopolymerPairs.get(varHomopolymerPairIndex);
            HomopolymerPair beforeHomopolymerPair = varHomopolymerPairIndex == 0 ? null : homopolymerPairs.get(varHomopolymerPairIndex - 1);
            HomopolymerPair afterHomopolymerPair = varHomopolymerPairIndex == homopolymerPairs.size() - 1 ? null : homopolymerPairs.get(varHomopolymerPairIndex + 1);

            if (varHomopolymerPair.base() != variant.Alt.charAt(0) || varHomopolymerPair.indelLength() != 1)
            {
                return false;
            }

            // check var pair and before pair
            if(beforeHomopolymerPair != null
                    && beforeHomopolymerPair.base() == variant.Ref.charAt(0)
                    && beforeHomopolymerPair.indelLength() == -1)
            {
                return true;
            }

            // check var pair and after pair
            if(afterHomopolymerPair != null
                    && afterHomopolymerPair.base() == variant.Ref.charAt(0)
                    && afterHomopolymerPair.indelLength() == -1)
            {
                return true;
            }

            return false;
        }
        else if(variant.indelLength() <= -1)
        {
            String deletedBases = variant.Ref.substring(1);
            List<Homopolymer> deletedAsHomopolymers = getHomopolymers(deletedBases);

            char deletedBase = deletedBases.charAt(0);
            int matchingPairCount = 0;
            for(int i = 0; i < homopolymerPairs.size(); ++i)
            {
                boolean matches = true;
                for(int j = 0; j < deletedAsHomopolymers.size(); ++j)
                {
                    Homopolymer deleted = deletedAsHomopolymers.get(j);
                    HomopolymerPair pair = homopolymerPairs.get(i + j);
                    if(deleted.Base != pair.base())
                    {
                        matches = false;
                        break;
                    }

                    if(-deleted.Length != pair.indelLength())
                    {
                        matches = false;
                        break;
                    }

                    if(j > 0 && j < deletedAsHomopolymers.size() - 1 && pair.ReadHomopolymer != null)
                    {
                        matches = false;
                        break;
                    }
                }

                if(matches)
                {
                    ++matchingPairCount;
                }
            }

            // TODO: not more than one match.
            assert matchingPairCount <= 1;

            return matchingPairCount == 1;
        }
        else
        {
            // TODO:
            throw new RuntimeException("TODO");
        }
    }

    private static int getVariantRefIndex(final VariantReadContext readContext)
    {
        // TODO: This is repetitive?
        String readCigar = readContext.readCigar();
        List<CigarOperator> cigarElements = cigarStringToOps(readCigar);
        List<CigarOperator> coreCigarOps = getCoreCigarOps(readContext, cigarElements);

        String coreRefBases = readContext.refBases();
        String coreReadBases = readContext.coreStr();
        int varCoreIndex = readContext.VarIndex - readContext.CoreIndexStart;

        int refIndex = 0;
        int readIndex = 0;
        for(CigarOperator op : coreCigarOps)
        {
            if(readIndex == varCoreIndex)
            {
                break;
            }

            if(op.consumesReferenceBases())
            {
                ++refIndex;
            }

            if(op.consumesReadBases())
            {
                ++readIndex;
            }
        }

        return refIndex;
    }

    private static boolean isSandwichedNonIndel(final VariantReadContext readContext, int varRefIndex)
    {
        if(readContext.variant().isIndel())
        {
            return false;
        }

        int varLength = readContext.variant().Ref.length();

        // TODO: Check that it isn't flanked by del or ins?
        // TODO: START this is repeat code
        String readCigar = readContext.readCigar();
        List<CigarOperator> cigarElements = cigarStringToOps(readCigar);
        List<CigarOperator> coreCigarOps = getCoreCigarOps(readContext, cigarElements);

        String coreRefBases = readContext.refBases();
        String coreReadBases = readContext.coreStr();
        int varCoreIndex = readContext.VarIndex - readContext.CoreIndexStart;

        int cigarIndex = 0;
        int readIndex = 0;
        for(CigarOperator op : coreCigarOps)
        {
            if(readIndex == varCoreIndex)
            {
                break;
            }

            if(op.consumesReadBases())
            {
                ++readIndex;
            }

            ++cigarIndex;
        }

        final List<CigarOperator> flankedVariantCigarOperators = coreCigarOps
                .subList(cigarIndex - 1, (cigarIndex - 1) + varLength + 2);

        assert flankedVariantCigarOperators.size() == varLength + 2;

        Set<CigarOperator> collapsedFlankedVariantCigarOps = flankedVariantCigarOperators
                .stream()
                .collect(Collectors.toCollection(() -> EnumSet.noneOf(CigarOperator.class)));

        assert collapsedFlankedVariantCigarOps.size() == 1;
        assert collapsedFlankedVariantCigarOps.stream().findFirst().orElse(null) == CigarOperator.M;
        // TODO: END

        String flankedRef = new String(readContext.RefBases, varRefIndex - 1, varLength + 2);
        String flankedRead = new String(readContext.ReadBases, readContext.VarIndex - 1, varLength + 2);

        assert flankedRef.charAt(0) == flankedRead.charAt(0);
        assert flankedRef.charAt(flankedRef.length() - 1) == flankedRead.charAt(flankedRead.length() - 1);

        if(flankedRef.charAt(0) != flankedRef.charAt(flankedRef.length() - 1))
        {
            return false;
        }

        char flankBase = flankedRef.charAt(0);
        String ref = readContext.variant().Ref;
        String alt = readContext.variant().Alt;
        String repeatedFlankBase = String.valueOf(flankBase).repeat(varLength);

        return ref.equals(repeatedFlankBase) || alt.equals(repeatedFlankBase);
    }

    // TODO: rename
    private static class AlignmentScore implements Comparable<AlignmentScore>
    {
        public final int FullIndelScore;
        public final int MatchScore;

        public AlignmentScore(final int fullIndelScore, final int matchScore)
        {
            FullIndelScore = fullIndelScore;
            MatchScore = matchScore;
        }

        public AlignmentScore incFullIndelScore(int i)
        {
            return new AlignmentScore(FullIndelScore + i, MatchScore);
        }

        public AlignmentScore incMatchScore(int i)
        {
            return new AlignmentScore(FullIndelScore, MatchScore + i);
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
            {
                return true;
            }
            if(!(o instanceof AlignmentScore))
            {
                return false;
            }
            final AlignmentScore that = (AlignmentScore) o;
            return FullIndelScore == that.FullIndelScore && MatchScore == that.MatchScore;
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(FullIndelScore, MatchScore);
        }

        @Override
        public String toString()
        {
            return format("(%d, %d)", FullIndelScore, MatchScore);
        }

        @Override
        public int compareTo(final AlignmentScore o)
        {
            int fullIndelScoreComparison = FullIndelScore - o.FullIndelScore;
            if(fullIndelScoreComparison != 0)
            {
                return fullIndelScoreComparison;
            }

            return MatchScore - o.MatchScore;
        }
    }

    private static enum AlignmentScoreTableDirection
    {
        UP(-1, 0),
        LEFT(0, -1),
        UP_LEFT(-1, -1);

        public final int RowInc;
        public final int ColInc;

        private AlignmentScoreTableDirection(int rowInc, int colInc)
        {
            RowInc = rowInc;
            ColInc = colInc;
        }
    }

    private static class AlignmentScoreTableElement
    {
        public final AlignmentScore Score;
        public final List<AlignmentScoreTableDirection> Directions;

        public AlignmentScoreTableElement(final AlignmentScore score, final List<AlignmentScoreTableDirection> directions)
        {
            Score = score;
            Directions = directions;
        }

        @Override
        public String toString()
        {
            return String.valueOf(Score) + ", " + Directions;
        }
    }

    private static void pairHomopolymersBacktrackHelper(final List<Homopolymer> refHomopolymers,
            final List<Homopolymer> readHomopolymers, final List<List<AlignmentScoreTableElement>> scores, int row, int col,
            final PVector<HomopolymerPair> currentAlignment, final List<List<HomopolymerPair>> alignmentsOut)
    {
        if(row == 0 && col == 0)
        {
            List<HomopolymerPair> finalAlignment = Lists.newArrayList();
            for(int i = currentAlignment.size() - 1; i >= 0; --i)
            {
                finalAlignment.add(currentAlignment.get(i));
            }

            alignmentsOut.add(finalAlignment);
            return;
        }

        if(row == 0)
        {
            PVector<HomopolymerPair> nextCurrentAlignment = currentAlignment.plus(new HomopolymerPair(null, readHomopolymers.get(col - 1)));
            pairHomopolymersBacktrackHelper(refHomopolymers, readHomopolymers, scores, row, col - 1, nextCurrentAlignment, alignmentsOut);
            return;
        }

        if(col == 0)
        {
            PVector<HomopolymerPair> nextCurrentAlignment = currentAlignment.plus(new HomopolymerPair(refHomopolymers.get(row - 1), null));
            pairHomopolymersBacktrackHelper(refHomopolymers, readHomopolymers, scores, row - 1, col, nextCurrentAlignment, alignmentsOut);
            return;
        }

        for(AlignmentScoreTableDirection direction : scores.get(row).get(col).Directions)
        {
            Homopolymer refHomopolymer = direction.RowInc == 0 ? null : refHomopolymers.get(row - 1);
            Homopolymer readHomopolymer = direction.ColInc == 0 ? null : readHomopolymers.get(col - 1);
            PVector<HomopolymerPair> nextCurrentAlignment;
            if(refHomopolymer == null || readHomopolymer == null || refHomopolymer.Base == readHomopolymer.Base)
            {
                nextCurrentAlignment = currentAlignment.plus(new HomopolymerPair(refHomopolymer, readHomopolymer));
            }
            else
            {
                nextCurrentAlignment = currentAlignment.plus(new HomopolymerPair(refHomopolymer, null)).plus(new HomopolymerPair(null, readHomopolymer));
            }

            pairHomopolymersBacktrackHelper(refHomopolymers, readHomopolymers, scores, row + direction.RowInc, col + direction.ColInc, nextCurrentAlignment, alignmentsOut);
        }
    }

    private static List<List<HomopolymerPair>> pairHomopolymersBacktrack(final List<Homopolymer> refHomopolymers,
            final List<Homopolymer> readHomopolymers, final List<List<AlignmentScoreTableElement>> scores)
    {
        int rowCount = refHomopolymers.size() + 1;
        int colCount = readHomopolymers.size() + 1;
        List<List<HomopolymerPair>> alignments = Lists.newArrayList();
        pairHomopolymersBacktrackHelper(refHomopolymers, readHomopolymers, scores, rowCount - 1, colCount - 1, TreePVector.empty(), alignments);
        return alignments;
    }

    // TODO: Do we need all of the pairings?
    @VisibleForTesting
    public static List<List<HomopolymerPair>> pairHomopolymers(final List<Homopolymer> refHomopolymers, final List<Homopolymer> readHomopolymers)
    {
        // we ignore the length of homopolymers and use the Needleman–Wunsch algorithm to pair up the ref and read homopolymers
        // we do not allow mismatches and simply penalise gaps
        // we are trying to minimize, not maximize.

        if(refHomopolymers.isEmpty() && readHomopolymers.isEmpty())
        {
            return Lists.newArrayList();
        }

        if(refHomopolymers.isEmpty())
        {
            return List.of(readHomopolymers.stream().map(homopolymer -> new HomopolymerPair(null, homopolymer)).collect(Collectors.toList()));
        }

        if(readHomopolymers.isEmpty())
        {
            return List.of(refHomopolymers.stream().map(homopolymer -> new HomopolymerPair(homopolymer,null)).collect(Collectors.toList()));
        }

        // now there is at least one homopolymer in both ref and read
        // construct scoring matrix, rows (resp. columns) correspond to ref (resp. read) homopolymers
        int rowCount = refHomopolymers.size() + 1;
        int colCount = readHomopolymers.size() + 1;
        List<List<AlignmentScoreTableElement>> scoreTableRows = Lists.newArrayList();
        scoreTableRows.add(IntStream.range(0, colCount).mapToObj(x -> new AlignmentScoreTableElement(new AlignmentScore(-x, 0), null)).collect(Collectors.toList()));

        List<AlignmentScore> tmpScores = Lists.newArrayList(null, null, null);
        for(int row = 1; row < rowCount; ++row)
        {
            Homopolymer refHomopolymer = refHomopolymers.get(row - 1);
            List<AlignmentScoreTableElement> nextRow = Lists.newArrayList(new AlignmentScoreTableElement(new AlignmentScore(-row, 0), null));
            for(int col = 1; col < colCount; ++col)
            {
                Homopolymer readHomopolymer = readHomopolymers.get(col - 1);

                // get the UP score, i.e. leaving off last ref homopolymer
                AlignmentScore upScore = scoreTableRows.get(row - 1).get(col).Score;
                tmpScores.set(AlignmentScoreTableDirection.UP.ordinal(), upScore.incFullIndelScore(-1));

                // get the LEFT score, i.e. leaving off last read homopolymer
                AlignmentScore leftScore = nextRow.get(col - 1).Score;
                tmpScores.set(AlignmentScoreTableDirection.LEFT.ordinal(), leftScore.incFullIndelScore(-1));

                // get the UP-LEFT score, i.e. attempting to match last ref and read homopolymer
                AlignmentScore upLeftScore = scoreTableRows.get(row - 1).get(col - 1).Score;
                if(refHomopolymer.Base == readHomopolymer.Base)
                {
                    int matchCount = min(refHomopolymer.Length, readHomopolymer.Length);
                    int mismatchCount = abs(refHomopolymer.Length - readHomopolymer.Length);
                    tmpScores.set(AlignmentScoreTableDirection.UP_LEFT.ordinal(), upLeftScore.incMatchScore(matchCount - mismatchCount));
                }
                else
                {
                    tmpScores.set(AlignmentScoreTableDirection.UP_LEFT.ordinal(), upLeftScore.incFullIndelScore(-2));
                }

                AlignmentScore bestScore = Collections.max(tmpScores);
                List<AlignmentScoreTableDirection> bestDirections = Lists.newArrayList();
                for(AlignmentScoreTableDirection direction : AlignmentScoreTableDirection.values())
                {
                    if(tmpScores.get(direction.ordinal()).equals(bestScore))
                    {
                        bestDirections.add(direction);
                    }
                }

                nextRow.add(new AlignmentScoreTableElement(bestScore, bestDirections));
            }

            scoreTableRows.add(nextRow);
        }

        // TODO: remove
        // Dump table.
        System.out.println("++++++++++++++++++++++++++++++");
        scoreTableRows.forEach(row -> System.out.println(row));
        System.out.println("------------------------------");

        return pairHomopolymersBacktrack(refHomopolymers, readHomopolymers, scoreTableRows);
    }

    private static List<CigarElement> collapseCigarOps(final List<CigarOperator> coreCigarOps)
    {
        List<CigarElement> cigarElements = Lists.newArrayList();
        CigarOperator currentOp = null;
        int length = 0;
        for(CigarOperator op : coreCigarOps)
        {
            if(currentOp == null)
            {
                currentOp = op;
                length = 1;
                continue;
            }

            if(currentOp == op)
            {
                ++length;
                continue;
            }

            cigarElements.add(new CigarElement(length, currentOp));
            currentOp = op;
            length = 1;
        }

        if(length >= 1)
        {
            cigarElements.add(new CigarElement(length, currentOp));
        }

        return cigarElements;
    }

    private static List<CigarOperator> getCoreCigarOps(final VariantReadContext readContext, final List<CigarOperator> cigarElements)
    {
        int cigarIndex = 0;
        int readIndex = 0;
        List<CigarOperator> coreCigarOps = Lists.newArrayList();
        while(readIndex <= readContext.CoreIndexEnd)
        {
            CigarOperator op = cigarElements.get(cigarIndex++);
            boolean inCore = readIndex >= readContext.CoreIndexStart;
            if(inCore)
            {
                coreCigarOps.add(op);
            }

            if(op.consumesReadBases())
            {
                ++readIndex;
            }
        }

        return coreCigarOps;
    }

    // TODO: Use library method.
    private static List<CigarOperator> cigarStringToOps(final String readCigar)
    {
        List<CigarOperator> cigarOps= Lists.newArrayList();
        int count = 0;
        for(int i = 0; i < readCigar.length(); ++i)
        {
            char c = readCigar.charAt(i);
            if('0' <= c && c <= '9')
            {
                count = 10*count + (int)c - (int)'0';
                continue;
            }

            if(count == 0)
            {
                throw new RuntimeException("TODO: Count should not be zero");
            }

            CigarOperator op;
            switch(c)
            {
                case 'M':
                    op = CigarOperator.M;
                    break;
                case 'D':
                    op = CigarOperator.D;
                    break;
                default:
                    throw new RuntimeException("TODO: Cigar operator not recognized");
            }

            for(; count > 0; --count)
            {
                cigarOps.add(op);
            }
        }

        if(count > 0)
        {
            throw new RuntimeException("TODO: Count should be zero");
        }

        return cigarOps;
    }
}
