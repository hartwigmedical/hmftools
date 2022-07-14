package com.hartwig.hmftools.common.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.DELIMITER;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.File;
import java.util.StringJoiner;

import org.jetbrains.annotations.NotNull;

public class AltSpliceJunctionFile
{
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

    public static final String FLD_ALT_SJ_POS_START = "SjStart";
    public static final String FLD_ALT_SJ_POS_END = "SjEnd";
    public static final String FLD_ALT_SJ_TYPE = "Type";
    public static final String FLD_ALT_SJ_FRAG_COUNT = "FragCount";

    public static final String ALT_SJ_FILE_ID = "alt_splice_junc.csv";

    public AltSpliceJunctionFile(
            final String geneId, final String geneName, final String chromosome, final int[] spliceJunction,
            final AltSpliceJunctionType type, final int fragmentCount, final int[] depthCounts,
            final AltSpliceJunctionContext[] regionContexts, final String[] baseContexts, final String[] transcriptNames)
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
    }

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + ISF_FILE_ID + ALT_SJ_FILE_ID;
    }

    public static String csvHeader()
    {
        return new StringJoiner(DELIMITER)
                .add(FLD_GENE_ID)
                .add(FLD_GENE_NAME)
                .add(FLD_CHROMOSOME)
                .add(FLD_ALT_SJ_POS_START)
                .add(FLD_ALT_SJ_POS_END)
                .add(FLD_ALT_SJ_TYPE)
                .add(FLD_ALT_SJ_FRAG_COUNT)
                .add("DepthStart")
                .add("DepthEnd")
                .add("RegionStart")
                .add("RegionEnd")
                .add("BasesStart")
                .add("BasesEnd")
                .add("TransStart")
                .add("TransEnd")
                .toString();
    }

    public boolean matches(final AltSpliceJunctionFile other)
    {
        return Chromosome.equals(other.Chromosome)
                && SpliceJunction[SE_START] == other.SpliceJunction[SE_START]
                && SpliceJunction[SE_END] == other.SpliceJunction[SE_END];
    }

    public int length() { return SpliceJunction[SE_END] - SpliceJunction[SE_START]; }

    public static String formKey(final String chromosome, int posStart, int posEnd)
    {
        return String.format("%s;%s;%s", chromosome, posStart, posEnd);
    }

    public String key() { return formKey(Chromosome, SpliceJunction[SE_START], SpliceJunction[SE_END]); }

    public String toCsv()
    {
        return new StringJoiner(DELIMITER)
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
                .add(TranscriptNames[SE_START])
                .add(TranscriptNames[SE_END])
                .toString();
    }

    public static AltSpliceJunctionFile fromCsv(final String data)
    {
        final String[] items = data.split(DELIMITER);

        int index = 0;

        String geneId = items[index++];
        String geneName = items[index++];
        String chromosome = items[index++];

        int[] spliceJunction = { Integer.parseInt(items[index++]), Integer.parseInt(items[index++]) };

        int fragCount = Integer.parseInt(items[index++]);

        int[] depthCounts = { Integer.parseInt(items[index++]), Integer.parseInt(items[index++]) };

        AltSpliceJunctionType type = AltSpliceJunctionType.valueOf(items[index++]);

        AltSpliceJunctionContext[] regionContexts =
                { AltSpliceJunctionContext.valueOf(items[index++]),
                        AltSpliceJunctionContext.valueOf(items[index++]) };

        String[] baseContexts = new String[] { items[index++], items[index++] };
        String[] transNames = new String[] { items[index++], items[index++] };

        return new AltSpliceJunctionFile(
                geneId, geneName, chromosome, spliceJunction, type,  fragCount, depthCounts,
                regionContexts, baseContexts, transNames);
    }

    public static AltSpliceJunctionFile fromCsv(
            final String[] items,
            int geneId, int geneName, int chr, int posStart, int posEnd, int type,
            int fragCount, int depthStart, int depthEnd, int regionStart, int regionEnd, int basesStart, int basesEnd, int transStart, int transEnd)
    {
        return new AltSpliceJunctionFile(
                items[geneId], items[geneName], items[chr], new int[] { Integer.parseInt(items[posStart]), Integer.parseInt(items[posEnd]) },
                AltSpliceJunctionType.valueOf(items[type]),
                Integer.parseInt(items[fragCount]), new int[] { Integer.parseInt(items[depthStart]), Integer.parseInt(items[depthEnd]) },
                new AltSpliceJunctionContext[] { AltSpliceJunctionContext.valueOf(items[regionStart]), AltSpliceJunctionContext.valueOf(items[regionEnd]) },
                new String[] { items[basesStart], items[basesEnd] }, new String[] { items[transStart], items[transEnd] } );
    }
}
