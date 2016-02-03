package com.hartwig.hmftools.sullivan;

import htsjdk.samtools.fastq.FastqRecord;
import org.jetbrains.annotations.NotNull;

public class SullivanRecord {

    @NotNull
    private final String header;
    @NotNull
    private final String sequence;
    @NotNull
    private final String qualityScores;

    @NotNull
    public static SullivanRecord createFromFastqRecord(@NotNull FastqRecord record, @NotNull FastqHeaderParser parser) {
        return new SullivanRecord(parser.apply(record.getReadHeader()),
                record.getReadString(), record.getBaseQualityString());
    }

    public SullivanRecord(@NotNull String header, @NotNull String sequence, @NotNull String qualityScores) {
        this.header = header;
        this.sequence = sequence;
        this.qualityScores = qualityScores;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SullivanRecord that = (SullivanRecord) o;

        return header.equals(that.header) && sequence.equals(that.sequence) && qualityScores.equals(that.qualityScores);
    }

    @Override
    public int hashCode() {
        int result = header.hashCode();
        result = 31 * result + sequence.hashCode();
        result = 31 * result + qualityScores.hashCode();
        return result;
    }

    @Override
    public String toString() {
        return "FastQ {" +
                "Header='" + header + '\'' +
                ", Sequence='" + sequence + '\'' +
                ", QualityScores='" + qualityScores + '\'' +
                '}';
    }
}
