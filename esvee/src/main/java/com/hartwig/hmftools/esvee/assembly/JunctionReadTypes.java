package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MAX_JUNC_POS_DIFF;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH_LOWER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.JUNCTION_PROXIMATE_READ_DISTANCE;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readIndexFromPosition;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.recordSoftClipsAtJunction;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.Junction;

public class JunctionReadTypes
{
    private final Junction mJunction;
    private final List<ReadJunctionInfo> mCandidateJunctionReads;
    private final List<Read> mNonJunctionReads;
    private final RefGenomeInterface mRefGenome;

    public JunctionReadTypes(final Junction junction, final List<Read> rawReads, final RefGenomeInterface refGenome)
    {
        mJunction = junction;
        mRefGenome = refGenome;

        mCandidateJunctionReads = Lists.newArrayList();
        mNonJunctionReads = Lists.newArrayList();

        for(Read read : rawReads)
        {
            ReadJunctionInfo readInfo = establishJunctionInfo(read);

            if(readInfo != null)
                mCandidateJunctionReads.add(readInfo);
            else
                mNonJunctionReads.add(read);
        }
    }

    public List<Read> nonJunctionReads() { return mNonJunctionReads; }

    public List<Read> junctionReads()
    {
        return mCandidateJunctionReads.stream().map(x -> x.Read).collect(Collectors.toList());
    }

    public void establishExtensionReads(final List<Read> extensionReads, final JunctionThresholdState juncState)
    {
        // apply standard and lower (SAGA-dependent) length limits in turn, or LINE if identified as a likely LINE assembly
        List<ReadJunctionInfo> candidates = mCandidateJunctionReads.stream().filter(x -> x.MatchesJunction).collect(Collectors.toList());
        int candidateCount = candidates.size();
        int lineCount = (int)candidates.stream().filter(x -> x.Read.hasLineTail(mJunction.isForward())).count();

        if(lineCount >= 0.5 * candidateCount)
        {
            juncState.IsLINE = true;
            juncState.MinExtensionLength = LINE_MIN_EXTENSION_LENGTH;
            juncState.MinSecondExtensionLength = LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH;
        }
        else
        {
            juncState.MinExtensionLength = ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
            juncState.MinSecondExtensionLength = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
        }

        for(ReadJunctionInfo readInfo : candidates)
        {
            if(readInfo.ExtensionLength < juncState.MinSecondExtensionLength)
                continue;

            if(!juncState.ExtensionLengthValid && readInfo.ExtensionLength >= juncState.MinExtensionLength)
            {
                juncState.ExtensionLengthValid = true;
            }
            else if(!juncState.SecondExtensionLengthValid)
            {
                juncState.SecondExtensionLengthValid = true;
            }

            extensionReads.add(readInfo.Read);
        }

        if(juncState.IsLINE)
            return;

        if(juncState.ExtensionLengthValid)
            return;

        // lower the limits and re-evaluate
        juncState.MinExtensionLength = ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER;
        juncState.MinSecondExtensionLength = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH_LOWER;
        juncState.UsesLowerSagaLimits = true;

        extensionReads.clear();

        for(ReadJunctionInfo readInfo : candidates)
        {
            if(readInfo.ExtensionLength < juncState.MinSecondExtensionLength)
                continue;

            if(!juncState.ExtensionLengthValid && readInfo.ExtensionLength >= juncState.MinExtensionLength)
            {
                juncState.ExtensionLengthValid = true;
            }
            else if(!juncState.SecondExtensionLengthValid)
            {
                juncState.SecondExtensionLengthValid = true;
            }

            extensionReads.add(readInfo.Read);
        }
    }

    private static final int MIN_SNV_TYPE_EXTENSION = 3; // ie to get at least 2 mismatches without clipping or indels

    private ReadJunctionInfo establishJunctionInfo(final Read read)
    {
        if(recordSoftClipsAtJunction(read, mJunction))
        {
            int extensionLength = readJunctionExtensionLength(read, mJunction);
            return new ReadJunctionInfo(read, extensionLength, true);
        }

        int extensionLength = 0;

        if(mJunction.isForward())
        {
            if(read.isRightClipped() && read.unclippedEnd() > mJunction.Position
            && abs(read.alignmentEnd() - mJunction.Position) <= ASSEMBLY_MAX_JUNC_POS_DIFF)
            {
                extensionLength = read.unclippedEnd() - mJunction.Position;
            }
            else if(read.hasIndelImpliedUnclippedEnd() && read.maxUnclippedEnd() > mJunction.Position) // indel-inferred soft-clips
            {
                extensionLength = read.maxUnclippedEnd() - mJunction.Position;
            }
            else if(!read.isRightClipped() && read.alignmentEnd() - mJunction.Position >= MIN_SNV_TYPE_EXTENSION
            && readMismatchesPastJunction(read, mJunction))
            {
                extensionLength = read.alignmentEnd() - mJunction.Position;
            }
        }
        else
        {
            if(read.isLeftClipped() && read.unclippedStart() < mJunction.Position
            && abs(read.alignmentStart() - mJunction.Position) <= ASSEMBLY_MAX_JUNC_POS_DIFF)
            {
                extensionLength = mJunction.Position - read.unclippedStart();
            }
            else if(read.hasIndelImpliedUnclippedStart() && read.minUnclippedStart() < mJunction.Position)
            {
                extensionLength = mJunction.Position - read.minUnclippedStart();
            }
            else if(!read.isLeftClipped() && mJunction.Position - read.alignmentStart() >= MIN_SNV_TYPE_EXTENSION
            && readMismatchesPastJunction(read, mJunction))
            {
                extensionLength = mJunction.Position - read.alignmentStart();
            }
        }

        return extensionLength > 0 ? new ReadJunctionInfo(read, extensionLength, false) : null;
    }

    @VisibleForTesting
    public boolean readMismatchesPastJunction(final Read read, final Junction junction)
    {
        int outerAlignmentPos = junction.isForward() ? read.alignmentEnd() : read.alignmentStart();

        int extensionLength = abs(outerAlignmentPos - junction.Position);
        if(extensionLength > ASSEMBLY_MAX_JUNC_POS_DIFF)
            return false;

        // must have a min number of mismatches past the junction
        Object numOfEvents = read.bamRecord().getAttribute(NUM_MUTATONS_ATTRIBUTE);

        if(numOfEvents == null || (int)numOfEvents < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
            return false;

        int posStart, posEnd;

        if(junction.isForward())
        {
            if(read.alignmentEnd() <= junction.Position)
                return false;

            posStart = junction.Position + 1;
            posEnd = read.alignmentEnd();
        }
        else
        {
            if(read.alignmentStart() >= junction.Position)
                return false;

            posStart = read.alignmentStart();
            posEnd = junction.Position - 1;
        }

        int refLength = posEnd - posStart + 1;

        if(refLength < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
            return false; // insufficient for junction support anyway

        if(refLength > MIN_VARIANT_LENGTH / 2)
            return false; // should have been soft-clipped

        byte[] refBases = mRefGenome.getBases(read.chromosome(), posStart, posEnd);

        if(refBases == null || refBases.length < refLength)
            return false;

        int readIndexStart = readIndexFromPosition(read, posStart);
        int mismatches = 0;

        for(int i = 0; i < refLength; ++i)
        {
            if(readIndexStart + i >= read.basesLength() || i >= refBases.length)
                return false;

            if(read.getBases()[readIndexStart + i] != refBases[i])
                ++mismatches;

            if(mismatches >= ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
                return true;
        }

        return false;
    }

    @VisibleForTesting
    public List<ReadJunctionInfo> candidateJunctionReads() { return mCandidateJunctionReads; }

    public static double calcProximateJunctionReadRatio(final Junction junction, final List<Read> rawReads)
    {
        if(!junction.softClipBased())
            return 0;

        // determine the ratio of extension reads vs any soft-clipped read within range - for use as a downstream filter
        int juncPosition = junction.Position;

        int proximateCount = 0;
        int exactExtensionCount = 0;

        for(Read read : rawReads)
        {
            if(junction.isForward())
            {
                if(read.isRightClipped() && read.unclippedEnd() > juncPosition
                && abs(read.alignmentEnd() - juncPosition) <= JUNCTION_PROXIMATE_READ_DISTANCE)
                {
                    ++proximateCount;

                    if(read.alignmentEnd() == juncPosition && read.unclippedEnd() - juncPosition >= ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
                        ++exactExtensionCount;
                }
            }
            else
            {
                if(read.isLeftClipped() && read.unclippedStart() < juncPosition
                && abs(read.alignmentStart() - juncPosition) <= JUNCTION_PROXIMATE_READ_DISTANCE)
                {
                    ++proximateCount;

                    if(read.alignmentStart() == juncPosition && juncPosition - read.unclippedStart() >= ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
                        ++exactExtensionCount;
                }
            }
        }

        return proximateCount > 0 ? exactExtensionCount / (double)proximateCount : 0;
    }
}
