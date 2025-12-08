package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;

import com.hartwig.hmftools.esvee.assembly.read.Read;

public class ReadAssemblyIndices
{
    public final int ReadIndexStart;
    public final int ReadIndexEnd;
    public final int AssemblyIndexStart;
    public final int JunctionIndex;

    public static final ReadAssemblyIndices INVALID_INDICES = new ReadAssemblyIndices(
            INVALID_INDEX, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);

    public ReadAssemblyIndices(final int readIndexStart, final int readIndexEnd, final int assemblyIndexStart)
    {
        this(readIndexStart, readIndexEnd, assemblyIndexStart, INVALID_INDEX);
    }

    public ReadAssemblyIndices(final int readIndexStart, final int readIndexEnd, final int assemblyIndexStart, final int junctionIndex)
    {
        ReadIndexStart = readIndexStart;
        ReadIndexEnd = readIndexEnd;
        AssemblyIndexStart = assemblyIndexStart;
        JunctionIndex = junctionIndex;
    }

    public boolean isValid() { return ReadIndexStart == INVALID_INDEX; }

    public int junctionReadStartDistance(int assemblyJunctionIndex)
    {
        // positive if the read's start is lower than then assembly's junction index, so will always be positive for junction reads
        // will only be negative for discordant reads and junction mates on -ve orientation assemblies
        return assemblyJunctionIndex - AssemblyIndexStart + ReadIndexStart;
    }

    public String toString()
    {
        return format("read(%d-%d) assembly(%d junc=%d)",
                ReadIndexStart, ReadIndexEnd, AssemblyIndexStart, JunctionIndex);
    }

    public static ReadAssemblyIndices getRefReadIndices(final JunctionAssembly assembly, final int refBasePosition, final Read read)
    {
        // ensure no indel-adjusted unclipped start is used when aligned to ref bases
        int junctionPosition = assembly.junction().Position;
        int alignmentStart = read.alignmentStart();
        int leftSoftClipOffset = read.leftClipLength();
        int rightSoftClipOffset = read.rightClipLength();
        int readIndexEnd = read.basesLength() - rightSoftClipOffset - 1;

        if(assembly.isForwardJunction())
        {
            if(alignmentStart >= junctionPosition) // not a ref read since past the junction
                return null;

            if(alignmentStart < refBasePosition)
            {
                // allow reads to start before the assembly's ref boundaries but will cull the bases used
                // TODO: is this to handle minor diffs due to soft-clips or if not, why didn't these reads extend the ref position?
                // suggest removing this allowance
                int readIndex = refBasePosition - alignmentStart + leftSoftClipOffset;
                return new ReadAssemblyIndices(readIndex, readIndexEnd, 0);
            }

            return new ReadAssemblyIndices(leftSoftClipOffset, readIndexEnd, alignmentStart - refBasePosition);
        }
        else
        {
            if(read.alignmentEnd() <= junctionPosition) // ends before the ref bases start
                return null;

            if(alignmentStart < junctionPosition)
            {
                // index off the relative start positions
                int readIndex = junctionPosition - alignmentStart + leftSoftClipOffset;
                return new ReadAssemblyIndices(readIndex, readIndexEnd, assembly.junctionIndex());
            }

            int assemblyIndex = assembly.junctionIndex() + alignmentStart - junctionPosition;
            return new ReadAssemblyIndices(leftSoftClipOffset, readIndexEnd, assemblyIndex);
        }
    }
}
