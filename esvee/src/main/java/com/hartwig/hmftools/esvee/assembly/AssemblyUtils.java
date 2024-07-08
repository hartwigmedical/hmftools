package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_DEL_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PROXIMATE_DUP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.NO_LINK;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.SECONDARY;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.SUPP_ONLY;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedBaseLength;
import static com.hartwig.hmftools.esvee.common.CommonUtils.createByteArray;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public final class AssemblyUtils
{
    public static int mismatchesPerComparisonLength(final int sequenceLength)
    {
        return (int)floor(log10(sequenceLength + 1));
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

    public static final byte N_BASE = 78;

    public static boolean basesMatch(
            final byte first, final byte second, final byte firstQual, final byte secondQual, final int lowQualThreshold)
    {
        return first == second || first == N_BASE || second == N_BASE
                || firstQual < lowQualThreshold || secondQual < lowQualThreshold;
    }

    public static boolean isLocalAssemblyCandidate(final JunctionAssembly first, final JunctionAssembly second)
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

        // must have concordant reads with mates crossing the other junction
        JunctionAssembly lowerAssembly = firstIsLower ? first : second;
        JunctionAssembly upperAssembly = !firstIsLower ? first : second;
        Junction lowerJunction = firstIsLower ? first.junction() : second.junction();
        Junction upperJunction = !firstIsLower ? first.junction() : second.junction();

        if(lowerAssembly.support().stream().noneMatch(x -> isCrossingConcordantRead(x, upperJunction, false)))
            return false;

        if(upperAssembly.support().stream().noneMatch(x -> isCrossingConcordantRead(x, lowerJunction, true)))
            return false;

        return true;
    }

    private static boolean isCrossingConcordantRead(final SupportRead read, final Junction junction, boolean requireLower)
    {
        if(read.isDiscordant() || read.isMateUnmapped() || !read.isPairedRead())
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

    public static byte[] createMinBaseQuals(final int length) { return createByteArray(length, (byte) (LOW_BASE_QUAL_THRESHOLD + 1)); }

    public static boolean hasUnsetBases(final JunctionAssembly assembly) { return !findUnsetBases(assembly.bases()).isEmpty(); }

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
}
