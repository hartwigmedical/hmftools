package com.hartwig.hmftools.common.genome.region;

public class ReplicationOriginRegion {

    public final String Chromosome;
    public final long Start;
    public final long End;
    public final double OriginValue;

    public ReplicationOriginRegion(final String chromosome, long start, long end, double originValue) {
        Chromosome = chromosome;
        Start = start;
        End = end;
        OriginValue = originValue;
    }
}
