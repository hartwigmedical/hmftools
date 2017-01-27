package com.hartwig.hmftools.retentionchecker;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.fastq.FastqRecord;

import org.jetbrains.annotations.NotNull;

class FastqHeader {

    private static final String PARSE_REGEXP = ":";

    @NotNull
    private final String instrumentName;
    private final int runId;
    @NotNull
    private final String flowcellId;
    private final int flowcellLane;
    @NotNull
    private final FastqHeaderKey key;

    @NotNull
    static FastqHeader parseFromFastqRecord(@NotNull final FastqRecord record,
            @NotNull final FastqHeaderNormalizer normalizer) {
        final String normalized = normalizer.apply(record.getReadHeader());

        final String[] parts = normalized.split(PARSE_REGEXP);
        return new FastqHeader(parts[0], Integer.valueOf(parts[1]), parts[2], Integer.valueOf(parts[3]),
                Integer.valueOf(parts[4]), Integer.valueOf(parts[5]), Integer.valueOf(parts[6]));
    }

    @VisibleForTesting
    FastqHeader(@NotNull final String instrumentName, final int runId, @NotNull final String flowcellId,
            final int flowcellLane, final int tileNumber, final int xCoordinate, final int yCoordinate) {
        this.instrumentName = instrumentName;
        this.runId = runId;
        this.flowcellId = flowcellId;
        this.flowcellLane = flowcellLane;
        this.key = new FastqHeaderKey(tileNumber, xCoordinate, yCoordinate);
    }

    @NotNull
    FastqHeaderKey key() {
        return key;
    }

    @NotNull
    String reference() {
        return instrumentName + PARSE_REGEXP + Integer.toString(runId) + PARSE_REGEXP + flowcellId + PARSE_REGEXP
                + Integer.toString(flowcellLane);
    }

    @SuppressWarnings("SimplifiableIfStatement")
    @Override
    public boolean equals(Object o) {
        if (this == o)
            return true;
        if (o == null || getClass() != o.getClass())
            return false;

        FastqHeader that = (FastqHeader) o;

        if (runId != that.runId)
            return false;
        if (flowcellLane != that.flowcellLane)
            return false;
        if (!instrumentName.equals(that.instrumentName))
            return false;
        if (!flowcellId.equals(that.flowcellId))
            return false;
        return key.equals(that.key);
    }

    @Override
    public int hashCode() {
        int result = instrumentName.hashCode();
        result = 31 * result + runId;
        result = 31 * result + flowcellId.hashCode();
        result = 31 * result + flowcellLane;
        result = 31 * result + key.hashCode();
        return result;
    }
}
