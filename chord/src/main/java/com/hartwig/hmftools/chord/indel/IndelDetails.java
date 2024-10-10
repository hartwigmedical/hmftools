package com.hartwig.hmftools.chord.indel;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

public class IndelDetails
{
    public final String mSampleId;

    public final String mChromosome;
    public final int mPosition;
    public final String mRefBases;
    public final String mAltBases;
    public final String mIndelSequence;
    public final int mIndelLength;
    public final String mRightFlankBases;
    public final String mLeftFlankBases;
    public final int mRightHomBasesCount;
    public final int mLeftHomBasesCount;
    public final int mMaxHomBasesCount;
    public final int mRepeatUnitsCount;

    public final MutationType mMutationType;
    public final ContextType mContextType;

    public static final int RIGHT_FLANK_INDEL_LENGTH_UNITS = 3;
    public static final int LEFT_FLANK_INDEL_LENGTH_UNITS = 1;

    private IndelDetails(
            String sampleId,
            String chromosome, int position, String refBases, String altBases,
            String indelSequence, int indelLength,
            String rightFlankBases, String leftFlankBases,
            int rightHomBasesCount, int leftHomBasesCount, int maxHomBasesCount, int repeatUnitsCount,
            MutationType mutationType, ContextType contextType
    )
    {
        mSampleId = sampleId;
        mChromosome = chromosome;
        mPosition = position;
        mRefBases = refBases;
        mAltBases = altBases;
        mIndelSequence = indelSequence;
        mIndelLength = indelLength;
        mRightFlankBases = rightFlankBases;
        mLeftFlankBases = leftFlankBases;
        mRightHomBasesCount = rightHomBasesCount;
        mLeftHomBasesCount = leftHomBasesCount;
        mMaxHomBasesCount = maxHomBasesCount;
        mRepeatUnitsCount = repeatUnitsCount;

        mMutationType = mutationType;
        mContextType = contextType;
    }

    public static IndelDetails from(String sampleId, IndelVariant indel, RefGenomeSource refGenome)
    {
        int indelLength = indel.mIndelLength;
        MutationType mutationType = indel.mIsDeletion ? MutationType.DELETION : MutationType.INSERTION;

        String indelSequence = indel.getIndelSequence();
        String rightFlankBases = indel.getFlankingRefBasesRight(refGenome, RIGHT_FLANK_INDEL_LENGTH_UNITS);
        String leftFlankBases = indel.getFlankingRefBasesLeft(refGenome, LEFT_FLANK_INDEL_LENGTH_UNITS);

        int rightHomBasesCount = IndelVariant.countHomologousBases(indelSequence, rightFlankBases);
        int leftHomBasesCount = IndelVariant.countHomologousBases(reverse(indelSequence), reverse(leftFlankBases));
        int maxHomBasesCount = Math.max(rightHomBasesCount, leftHomBasesCount);

        // Only need to check repeat units from the 5' -> 3' direction, since the reported POS of an indel in a repeat region is
        // always the first position of the repeat region.
        int repeatUnitsCount = IndelVariant.countRepeatUnits(indelSequence, rightFlankBases);

        ContextType contextType;
        if(repeatUnitsCount >= 2)
        {
            contextType = (indelLength < 50) ?
                    ContextType.REPEAT :
                    ContextType.MICROHOMOLOGY;
        }
        else if(repeatUnitsCount >= 1 && maxHomBasesCount >= 2)
        { // Require more homologous bases for shorter indels
            contextType = ContextType.MICROHOMOLOGY;
        }
        else if(repeatUnitsCount >= 1 && maxHomBasesCount >= 1 && indelLength > 3)
        { // Tolerate fewer homologous bases for longer indels
            contextType = ContextType.MICROHOMOLOGY;
        }
        else
        {
            contextType = ContextType.NONE;
        }

        return new IndelDetails(
                sampleId,
                indel.mVariant.Chromosome.toString(),
                indel.mVariant.Position,
                indel.mVariant.RefBases,
                indel.mVariant.AltBases,
                indelSequence, indelLength,
                rightFlankBases, leftFlankBases,
                rightHomBasesCount, leftHomBasesCount, maxHomBasesCount, repeatUnitsCount,
                mutationType, contextType
        );
    }

    private static String reverse(String sequence) { return new StringBuilder(sequence).reverse().toString(); }

    public static BufferedWriter initializeWriter(String path) throws IOException
    {
        BufferedWriter writer = FileWriterUtils.createBufferedWriter(path, false);

        String header = String.join(
                TSV_DELIM,
                "SampleId",
                "Chromosome",
                "Position",
                "RefBases",
                "AltBases",
                "IndelSequence",
                "IndelLength",
                "RightFlankBases",
                "LeftFlankBases",
                "RightHomBasesCount",
                "LeftHomBasesCount",
                "MaxHomBasesCount",
                "RepeatUnitsCount",
                "MutationType",
                "ContextType"
        );

        writer.write(header);
        writer.newLine();

        return writer;
    }

    public void writeLine(BufferedWriter writer) throws IOException
    {
        String line = String.join(
                TSV_DELIM,
                mSampleId,
                mChromosome,
                String.valueOf(mPosition),
                mRefBases,
                mAltBases,
                mIndelSequence,
                String.valueOf(mIndelLength),
                mRightFlankBases,
                mLeftFlankBases,
                String.valueOf(mRightHomBasesCount),
                String.valueOf(mLeftHomBasesCount),
                String.valueOf(mMaxHomBasesCount),
                String.valueOf(mRepeatUnitsCount),
                String.valueOf(mMutationType),
                String.valueOf(mContextType)
        );

        writer.write(line);
        writer.newLine();
    }
}
