package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsOverlap;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.sv.SvRegion;

public class LineInsertSiteData
{
    public final String Program;
    public final String SampleId;

    public final SvRegion InsertSite;
    public final String InsertType;
    public final SvRegion SourceSite;
    public final boolean HasInversion;

    // Linx only
    public final String ChainDesc;
    public final int ClusterId;
    public final int ChainId;

    public static final String INSERT_TYPE_SOLO_L1 = "SOLO_L1";
    public static final String INSERT_TYPE_PARTNERED = "PARTNERED";
    public static final String INSERT_TYPE_TRANSDUCTION = "TRANSDUCTION";

    public static final String PROGRAM_LINX = "LINX";
    public static final String PROGRAM_PCAWG = "PCAWG";

    public LineInsertSiteData(final String program, final String sampleId, final SvRegion insertSite, final String insertType,
            final SvRegion sourceSite, final boolean hasInversion, final int clusterId, final int chainId, final String chainDesc)
    {
        Program = program;
        SampleId = sampleId;
        InsertSite = insertSite;
        InsertType = insertType;
        SourceSite = sourceSite;
        HasInversion = hasInversion;
        ChainDesc = chainDesc;
        ClusterId = clusterId;
        ChainId = chainId;
    }

    private static final int MAX_BASE_DIFF = 50;

    public boolean matches(final LineInsertSiteData other)
    {
        if(!other.InsertSite.Chromosome.equals(InsertSite.Chromosome))
            return false;

        if(other.InsertSite.end() > 0 && InsertSite.end() > 0)
        {
            return positionsOverlap(
                    other.InsertSite.start(), other.InsertSite.end(),
                    InsertSite.start() - MAX_BASE_DIFF, InsertSite.end() + MAX_BASE_DIFF);
        }
        else
        {
            return abs(other.InsertSite.start() - InsertSite.start()) < MAX_BASE_DIFF;
        }
    }
}
