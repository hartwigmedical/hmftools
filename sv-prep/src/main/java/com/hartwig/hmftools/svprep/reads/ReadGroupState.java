package com.hartwig.hmftools.svprep.reads;

import static java.lang.String.format;

public class ReadGroupState
{
    public final String ReadId;
    public final String SourceChrPartition;
    public final String RemoteChrPartition;
    public final int RemotePosition;
    public final int ReadCount;
    public final ReadGroupStatus Status;
    public final boolean FirstHasSupplementary;
    public final boolean SecondHasSupplementary;

    public ReadGroupState(
            final String readId, final String sourceChrPartition, final String remoteChrPartition, final int remotePosition,
            final int readCount, final ReadGroupStatus status, final boolean firstHasSupplementary, final boolean secondHasSupplementary)
    {
        ReadId = readId;
        SourceChrPartition = sourceChrPartition;
        RemoteChrPartition = remoteChrPartition;
        RemotePosition = remotePosition;
        ReadCount = readCount;
        Status = status;
        FirstHasSupplementary = firstHasSupplementary;
        SecondHasSupplementary = secondHasSupplementary;
    }

    public String toString()
    {
        return format("id(%s) reads(%d) chrPartions(%s - %s) status(%s)", ReadId, ReadCount, SourceChrPartition, RemoteChrPartition, Status);
    }
}
