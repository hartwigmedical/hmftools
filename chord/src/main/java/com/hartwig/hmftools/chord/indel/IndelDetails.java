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

    // Thresholds for determining the context of an indel
    private static final int REPETITIVE_REGION_MIN_REPEAT_UNITS = 2;
    private static final int SHORT_REPEAT_UNIT_MAX_LENGTH = 50;

    private static final int MIN_HOMOLOGOUS_BASES = 2;
    private static final int LONG_INDEL_MIN_HOMOLOGOUS_BASES = 1;
    private static final int LONG_INDEL_MIN_LENGTH = 4;

    // Left flank only used to determine homology, therefore do not need to scan far
    private static final int LEFT_FLANK_INDEL_LENGTH_UNITS = 1;

    // Right flank used to determine homology, but also the number of repeat lengths after which we consider the indel region repetitive.
    // We therefore only need to scan just further than the number of repeat lengths after which we consider the indel region repetitive.
    private static final int RIGHT_FLANK_INDEL_LENGTH_UNITS = REPETITIVE_REGION_MIN_REPEAT_UNITS + 1; //

    public static IndelDetails from(String sampleId, IndelVariant indel, RefGenomeSource refGenome)
    {
        int indelLength = indel.mIndelLength;
        MutationType mutationType = indel.mIsDeletion ? MutationType.DELETION : MutationType.INSERTION;

        String indelSequence = indel.getIndelSequence();
        String rightFlankBases = indel.getFlankingRefBasesRight(refGenome, RIGHT_FLANK_INDEL_LENGTH_UNITS);
        String leftFlankBases = indel.getFlankingRefBasesLeft(refGenome, LEFT_FLANK_INDEL_LENGTH_UNITS);

        int rightHomBasesCount = IndelVariant.countHomologousBases(indelSequence, rightFlankBases);
        int leftHomBasesCount = IndelVariant.countHomologousBases(reverse(indelSequence), reverse(leftFlankBases));
        int homBasesCount = Math.max(rightHomBasesCount, leftHomBasesCount);

        // Only need to check repeat units from the 5' -> 3' direction, since the reported POS of an indel in a repeat region is
        // always the first position of the repeat region.
        int repeatUnitsCount = IndelVariant.countRepeatUnits(indelSequence, rightFlankBases);

        ContextType contextType = determineIndelContext(indelLength, repeatUnitsCount, homBasesCount);

        return new IndelDetails(
                sampleId,
                indel.mVariant.Chromosome,
                indel.mVariant.Position,
                indel.mVariant.RefBases,
                indel.mVariant.AltBases,
                indelSequence, indelLength,
                rightFlankBases, leftFlankBases,
                rightHomBasesCount, leftHomBasesCount, homBasesCount, repeatUnitsCount,
                mutationType, contextType
        );
    }

    private static ContextType determineIndelContext(int indelLength, int repeatUnitsCount, int homBasesCount)
    {
        ContextType contextType;

        if(repeatUnitsCount >= REPETITIVE_REGION_MIN_REPEAT_UNITS)
        {
            if(indelLength < SHORT_REPEAT_UNIT_MAX_LENGTH)
            {
                // In repetitive regions where the repeat unit is a short sequence, an indel is more likely due to polymerase slippage
                // rather than homology directed repair
                return ContextType.REPEAT;
            }
            else
            {
                // An indel with a long sequence exactly matching the flanking region is more likely a result of homology directed repair
                contextType = ContextType.MICROHOMOLOGY;
            }
        }
        else if(homBasesCount >= MIN_HOMOLOGOUS_BASES)
        {
            contextType = ContextType.MICROHOMOLOGY;
        }
        else if(homBasesCount >= LONG_INDEL_MIN_HOMOLOGOUS_BASES && indelLength >= LONG_INDEL_MIN_LENGTH)
        { // Tolerate fewer homologous bases for longer indels
            contextType = ContextType.MICROHOMOLOGY;
        }
        else
        {
            contextType = ContextType.NONE;
        }

        return contextType;
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
