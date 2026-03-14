package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;

public class AltSpliceJuncData
{
    // has more fields than the commonly used NovelSpliceJunction, used internally by Isofox for cohort and other analyses
    public final String GeneId;
    public final String GeneName;
    public final String Chromosome;
    public final int[] SpliceJunction;
    public final AltSpliceJunctionType Type;

    public final AltSpliceJunctionContext[] RegionContexts;

    public final int FragmentCount;
    public final int[] DepthCounts; // counts at the start and end
    public final String[] TranscriptNames;
    public final String[] BaseContexts;
    public final int CohortFrequency;

    public AltSpliceJuncData(
            final String geneId, final String geneName, final String chromosome, final int[] spliceJunction,
            final AltSpliceJunctionType type, final int fragmentCount, final int[] depthCounts,
            final AltSpliceJunctionContext[] regionContexts, final String[] baseContexts, final String[] transcriptNames,
            final int cohortFrequency)
    {
        GeneId = geneId;
        GeneName = geneName;
        Chromosome = chromosome;
        SpliceJunction = spliceJunction;
        Type = type;
        RegionContexts = regionContexts;
        FragmentCount = fragmentCount;
        DepthCounts = depthCounts;
        TranscriptNames = transcriptNames;
        BaseContexts = baseContexts;
        CohortFrequency = cohortFrequency;
    }

    public boolean matches(final AltSpliceJuncData other)
    {
        return Chromosome.equals(other.Chromosome)
                && SpliceJunction[SE_START] == other.SpliceJunction[SE_START]
                && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public int length() { return SpliceJunction[SE_END] - SpliceJunction[SE_START]; }

    public String key() { return formKey(Chromosome, SpliceJunction[SE_START], SpliceJunction[SE_END]); }
}
