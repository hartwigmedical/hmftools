package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_DEL_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PROXIMATE_DUP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.NO_LINK;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.SECONDARY;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.SUPP_ONLY;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedBaseLength;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.createByteArray;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public final class AssemblyUtils
{
    public static final int DNA_BASE_COUNT = Nucleotides.DNA_BASES.length + 1; // allows for Ns
    public static final byte NO_BASE = 0;

    public static int baseIndex(final byte base)
    {
        // protects against out of array errors from non-standard letters (N is permitted)
        int baseIndex = Nucleotides.baseIndex(base);
        return baseIndex < 0 || baseIndex >= DNA_BASE_COUNT ? DNA_BASE_COUNT - 1 : baseIndex;
    }

    public static int mismatchesPerComparisonLength(final int sequenceLength)
    {
        if(sequenceLength < 15)
            return 0;

        if(sequenceLength <= 100)
            return 1;

        return (int)ceil(sequenceLength / 200.0) + 1;
    }

    public static int readQualFromJunction(final Read read, final Junction junction)
    {
        int readJunctionIndex = read.getReadIndexAtReferencePosition(junction.Position, true);

        int readIndexStart;
        int readIndexEnd;

        if(junction.isForward())
        {
            readIndexStart = readJunctionIndex;
            readIndexEnd = read.basesLength() - 1;
        }
        else
        {
            readIndexStart = 0;
            readIndexEnd = readJunctionIndex;
        }

        readIndexStart = max(readIndexStart, 0);
        readIndexEnd = min(readIndexEnd, read.getBaseQuality().length - 1);

        int baseQualTotal = 0;

        for(int i = readIndexStart; i <= readIndexEnd; ++i)
        {
            baseQualTotal += read.getBaseQuality()[i];
        }

        return baseQualTotal;
    }

    public static boolean basesMatch(
            final byte first, final byte second, final byte firstQual, final byte secondQual)
    {
        return first == second || first == DNA_N_BYTE || second == DNA_N_BYTE || belowMinQual(firstQual) || belowMinQual(secondQual);
    }

    public static boolean basesMatch(final byte[] bases1, final byte[] bases2)
    {
        if(bases1.length != bases2.length)
            return false;

        for(int i = 0; i < bases1.length; ++i)
        {
            if(bases1[i] != bases2[i])
                return false;
        }

        return true;
    }

    public static boolean isLocalAssemblyCandidate(
            final JunctionAssembly first, final JunctionAssembly second, boolean checkConcordantReads, boolean checkLineInsertion)
    {
        if(!first.junction().Chromosome.equals(second.junction().Chromosome))
            return false;

        // assemblies must have DEL or DUP orientations, be within threshold distances of each other
        if(first.isForwardJunction() == second.isForwardJunction())
            return false;

        boolean firstIsLower = first.junction().Position <= second.junction().Position;
        boolean isDelType = firstIsLower == first.isForwardJunction();
        int junctionDistance = abs(first.junction().Position - second.junction().Position);

        if((isDelType && junctionDistance > PROXIMATE_DEL_LENGTH) || (!isDelType && junctionDistance > PROXIMATE_DUP_LENGTH))
            return false;

        if(!checkConcordantReads && !checkLineInsertion)
            return true;

        if(checkLineInsertion && (first.hasLineSequence() || second.hasLineSequence()))
            return true;

        if(checkConcordantReads)
        {
            // must have concordant reads with mates crossing the other junction
            JunctionAssembly lowerAssembly = firstIsLower ? first : second;
            JunctionAssembly upperAssembly = !firstIsLower ? first : second;
            Junction lowerJunction = firstIsLower ? first.junction() : second.junction();
            Junction upperJunction = !firstIsLower ? first.junction() : second.junction();

            if(lowerAssembly.support().stream().anyMatch(x -> isCrossingConcordantRead(x, upperJunction, false))
            || upperAssembly.support().stream().anyMatch(x -> isCrossingConcordantRead(x, lowerJunction, true)))
            {
                return true;
            }
        }

        return false;
    }

    private static boolean isCrossingConcordantRead(final SupportRead read, final Junction junction, boolean requireLower)
    {
        if(read.isDiscordant() || !read.isPairedRead() || read.isMateUnmapped())
            return false;

        if(requireLower)
            return read.mateAlignmentEnd() <= junction.Position;
        else
            return read.mateAlignmentStart() >= junction.Position;
    }

    public static boolean isSupplementaryOnly(final JunctionAssembly assembly)
    {
        return assembly.support().stream().allMatch(x -> x.isSupplementary());
    }

    public static byte[] createBaseQualsAboveMinThreshold(final int length) { return createByteArray(length, LOW_BASE_QUAL_THRESHOLD); }

    public static String extractInsertSequence(
            final JunctionAssembly first, boolean firstReversed, final JunctionAssembly second, boolean secondReversed, int insertLength)
    {
        int extBaseIndexStart, extBaseIndexEnd;

        if(first.isForwardJunction())
        {
            extBaseIndexStart = first.junctionIndex() + 1;
            extBaseIndexEnd = min(extBaseIndexStart + insertLength - 1, first.baseLength() - 1);
        }
        else
        {
            extBaseIndexEnd = first.junctionIndex() - 1;
            extBaseIndexStart = max(extBaseIndexEnd - insertLength + 1, 0);
        }

        String insertSequence = first.formSequence(extBaseIndexStart, extBaseIndexEnd);

        if(firstReversed)
            insertSequence = Nucleotides.reverseComplementBases(insertSequence);

        if(insertLength <= first.extensionLength())
            return insertSequence;

        // take the extra bases from the second assembly
        int remainingInsertLength = insertLength - first.extensionLength();

        if(second.isForwardJunction())
        {
            extBaseIndexStart = second.junctionIndex() + 1;
            extBaseIndexEnd = min(extBaseIndexStart + remainingInsertLength - 1, second.baseLength() - 1);
        }
        else
        {
            extBaseIndexEnd = second.junctionIndex() - 1;
            extBaseIndexStart = max(extBaseIndexEnd - remainingInsertLength + 1, 0);
        }

        String extraInsertSequence = second.formSequence(extBaseIndexStart, extBaseIndexEnd);

        if(secondReversed)
            extraInsertSequence = Nucleotides.reverseComplementBases(extraInsertSequence);

        return insertSequence + extraInsertSequence;
    }

    public static int calcTrimmedRefBaseLength(final JunctionAssembly assembly)
    {
        int refBaseLength = assembly.refBaseLength();

        if(assembly.repeatInfo().isEmpty())
            return refBaseLength;

        int seqStart = assembly.isForwardJunction() ? assembly.junctionIndex() - refBaseLength + 1 : assembly.junctionIndex();
        int seqEnd = assembly.isForwardJunction() ? assembly.junctionIndex() : assembly.junctionIndex() + refBaseLength - 1;

        return calcTrimmedBaseLength(seqStart, seqEnd, assembly.repeatInfo());
    }

    public static int calcTrimmedExtensionBaseLength(final JunctionAssembly assembly)
    {
        int extBaseLength = assembly.extensionLength();

        if(assembly.repeatInfo().isEmpty())
            return extBaseLength;

        int seqStart = assembly.isForwardJunction() ? assembly.junctionIndex() + 1 : 0;
        int seqEnd = assembly.isForwardJunction() ? seqStart + extBaseLength - 1 : extBaseLength - 1;

        return calcTrimmedBaseLength(seqStart, seqEnd, assembly.repeatInfo());
    }

    public static JunctionAssembly findMatchingAssembly(
            final List<JunctionAssembly> assemblies, final JunctionAssembly assembly, boolean requireExtensionMatch)
    {
        return assemblies.stream()
                .filter(x -> x != assembly)
                .filter(x -> x.junction().compareTo(assembly.junction()) == 0)
                .filter(x -> !requireExtensionMatch || x.extensionLength() == assembly.extensionLength())
                .findFirst().orElse(null);
    }

    public static void setAssemblyOutcome(final JunctionAssembly assembly)
    {
        if(assembly.outcome() != AssemblyOutcome.UNSET)
            return;

        if(assembly.phaseGroup() == null)
        {
            assembly.setOutcome(NO_LINK);
            return;
        }

        List<AssemblyLink> secondarySplitLinks = assembly.phaseGroup().findSecondarySplitLinks(assembly);

        if(!secondarySplitLinks.isEmpty())
        {
            assembly.setOutcome(SECONDARY);
            return;
        }

        // check for assemblies only comprised of supp reads for support
        if(isSupplementaryOnly(assembly))
        {
            assembly.setOutcome(SUPP_ONLY);
            return;
        }

        assembly.setOutcome(NO_LINK);
    }

    public static boolean hasUnsetBases(final JunctionAssembly assembly) { return !findUnsetBases(assembly.bases()).isEmpty(); }

    public static List<int[]> findUnsetBases(final byte[] bases)
    {
        List<int[]> emptyRanges = Lists.newArrayList();

        int[] range = null;

        for(int i = 0; i < bases.length; ++i)
        {
            if(bases[i] == 0)
            {
                if(range == null)
                {
                    range = new int[] {i, -1};
                    emptyRanges.add(range);
                }
            }
            else
            {
                if(range != null)
                {
                    range[1] = i - 1;
                    range = null;
                }
            }
        }

        if(range != null)
            range[1] = bases.length - 1;

        return emptyRanges;
    }

    public static String nonNullBaseStr(final byte[] bases)
    {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < bases.length; ++i)
        {
            if(bases[i] != 0)
                sb.append((char)bases[i]);
        }

        return sb.toString();
    }
}
