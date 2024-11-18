package com.hartwig.hmftools.chord.snv;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.chord.prep.SmallVariant;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.sigs.SnvSigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

public class SnvDetails
{
    public final String mSampleId;

    public final String mChromosome;
    public final int mPosition;
    public final String mRefBases;
    public final String mAltBases;

    public final String mTriNucSequence;
    public final String mTriNucContext;
    public final String mTriNucContextRenamed;

    private SnvDetails(
            String sampleId,
            String chromosome, int position, String refBases, String altBases,
            String triNucSequence, String triNucContext, String triNucContextRenamed
    )
    {
        mSampleId = sampleId;

        mChromosome = chromosome;
        mPosition = position;
        mRefBases = refBases;
        mAltBases = altBases;

        mTriNucSequence = triNucSequence;
        mTriNucContext = triNucContext;
        mTriNucContextRenamed = triNucContextRenamed;
    }

    public static SnvDetails from(String sampleId, SmallVariant snv, RefGenomeSource refGenome)
    {
        String triNucSequence = new String(refGenome.getBases(snv.Chromosome, snv.Position-1,  snv.Position+1));
        String triNucContext = SnvSigUtils.variantContext(snv.RefBases, snv.AltBases, triNucSequence);

        return new SnvDetails(
                sampleId,
                snv.Context.getContig(),
                snv.Context.getStart(),
                snv.RefBases,
                snv.AltBases,
                triNucSequence,
                triNucContext,
                renameTriNucBin(triNucContext)
        );
    }

    public static String renameTriNucBin(String triNucContext)
    {
        // Convert e.g. "C>A_ACA" to "A[C>A]A"
        // The latter is the format required by the CHORD random forest model

        String[] bucketNameSplit = triNucContext.split("_");

        String substitutionType = bucketNameSplit[0];
        char[] triNucSequence = bucketNameSplit[1].toCharArray();

        return String.format("%s[%s]%s", triNucSequence[0], substitutionType, triNucSequence[2]);
    }

    public static BufferedWriter initializeWriter(String path) throws IOException
    {
        BufferedWriter writer = FileWriterUtils.createBufferedWriter(path, false);

        String header = String.join(
                TSV_DELIM,
                "Chromosome",
                "Position",
                "RefBases",
                "AltBases",
                "TriNucSequence",
                "TriNucContext",
                "TriNucContextRenamed"
        );

        writer.write(header);
        writer.newLine();

        return writer;
    }

    public void writeLine(BufferedWriter writer) throws IOException
    {
        String line = String.join(
                TSV_DELIM,
                mChromosome,
                String.valueOf(mPosition),
                mRefBases,
                mAltBases,
                mTriNucSequence,
                mTriNucContext,
                mTriNucContextRenamed
        );

        writer.write(line);
        writer.newLine();
    }
}
