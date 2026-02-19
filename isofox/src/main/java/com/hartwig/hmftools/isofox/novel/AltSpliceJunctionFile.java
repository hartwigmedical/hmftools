package com.hartwig.hmftools.isofox.novel;

import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;

public class AltSpliceJunctionFile
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

    public static final String FLD_TRANS_START = "TransStart";
    public static final String FLD_TRANS_END = "TransEnd";


    public AltSpliceJunctionFile(
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

    public static String header()
    {
        // expands upon the fields required to load NovelSpliceJunction
        StringJoiner header = new StringJoiner(TSV_DELIM);
        header.add(FLD_GENE_ID);
        header.add(NovelSpliceJunctionFile.header());
        header.add(FLD_TRANS_START);
        header.add(FLD_TRANS_END);

        return header.toString();
    }

    public boolean matches(final AltSpliceJunctionFile other)
    {
        return Chromosome.equals(other.Chromosome)
                && SpliceJunction[SE_START] == other.SpliceJunction[SE_START]
                && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public int length() { return SpliceJunction[SE_END] - SpliceJunction[SE_START]; }

    public String key() { return formKey(Chromosome, SpliceJunction[SE_START], SpliceJunction[SE_END]); }

    public String toLine()
    {
        return new StringJoiner(TSV_DELIM)
                .add(GeneId)
                .add(GeneName)
                .add(Chromosome)
                .add(String.valueOf(SpliceJunction[SE_START]))
                .add(String.valueOf(SpliceJunction[SE_END]))
                .add(String.valueOf(Type))
                .add(String.valueOf(FragmentCount))
                .add(String.valueOf(DepthCounts[SE_START]))
                .add(String.valueOf(DepthCounts[SE_END]))
                .add(String.valueOf(RegionContexts[SE_START]))
                .add(String.valueOf(RegionContexts[SE_END]))
                .add(BaseContexts[SE_START])
                .add(BaseContexts[SE_END])
                .add(String.valueOf(CohortFrequency))
                .add(TranscriptNames[SE_START])
                .add(TranscriptNames[SE_END])
                .toString();
    }
}
