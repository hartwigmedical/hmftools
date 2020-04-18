package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.SPLICED_BOTH;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.SPLICED_ONE;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.UNKNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.lowerChromosome;

import java.util.List;

import com.hartwig.hmftools.isofox.common.ReadRecord;

import htsjdk.samtools.CigarOperator;

public class FusionFragment
{
    private final List<ReadRecord> mReads;

    private final String[] mChromosomes;
    private final long[] mSjPositions;
    private final byte[] mSjOrientations;
    private final boolean[] mSjValid;
    private FusionFragmentType mType;

    public FusionFragment(final List<ReadRecord> reads)
    {
        mReads = reads;

        mSjPositions = new long[] {-1, -1};
        mChromosomes = new String[]{"", ""};
        mSjOrientations = new byte[]{0, 0};
        mSjValid = new boolean[]{false, false};
        int sjCount = 0;

        for(final ReadRecord read : reads)
        {
            if(!read.Cigar.containsOperator(CigarOperator.S))
                continue;

            final List<long[]> readCoords = read.getMappedRegionCoords();

            boolean useLeft;
            long sjPosition = 0;
            byte sjOrientation = 0;

            if(read.Cigar.isLeftClipped() && read.Cigar.isRightClipped())
            {
                // should be very unlikely since implies a very short exon and even then would expect it to be mapped
                useLeft = read.Cigar.getFirstCigarElement().getLength() > read.Cigar.getLastCigarElement().getLength();
            }
            else
            {
                useLeft = read.Cigar.isLeftClipped();
            }

            if(useLeft)
            {
                sjPosition = readCoords.get(0)[SE_START];
                sjOrientation = -1;
            }
            else
            {
                sjPosition = readCoords.get(readCoords.size() - 1)[SE_END];
                sjOrientation = 1;
            }

            if(sjCount == 0)
            {
                mChromosomes[SE_START] = read.Chromosome;
                mSjPositions[SE_START] = sjPosition;
                mSjOrientations[SE_START] = sjOrientation;
                mSjValid[SE_START] = true;
            }
            else
            {
                if((mChromosomes[SE_START].equals(read.Chromosome) && mSjPositions[SE_START] < sjPosition)
                || (!mChromosomes[SE_START].equals(read.Chromosome) && lowerChromosome(mChromosomes[SE_START], read.Chromosome)))
                {
                    // already in correct positions
                    mChromosomes[SE_END] = read.Chromosome;
                    mSjPositions[SE_END] = sjPosition;
                    mSjOrientations[SE_END] = sjOrientation;
                    mSjValid[SE_END] = true;
                }
                else
                {
                    mChromosomes[SE_END] = mChromosomes[SE_START];
                    mSjPositions[SE_END] = mSjPositions[SE_START];
                    mSjOrientations[SE_END] = mSjOrientations[SE_START];
                    mSjValid[SE_END] = mSjValid[SE_START];

                    mChromosomes[SE_START] = read.Chromosome;
                    mSjPositions[SE_START] = sjPosition;
                    mSjOrientations[SE_START] = sjOrientation;
                    mSjValid[SE_START] = true;
                }
            }

            break;
        }

        mType = UNKNOWN;

            if(mSjValid[SE_START] && mSjValid[SE_END])
        {
            mType = SPLICED_BOTH;
        }
        else if(mSjValid[SE_START] || mSjValid[SE_END])
        {
            mType = SPLICED_ONE;
        }
        else
        {
            mType = DISCORDANT;
        }
    }

    public FusionFragmentType type() { return mType; }
    public final long[] splicePositions() { return mSjPositions; }
    public final byte[] spliceOrientations() { return mSjOrientations; }
    public final String[] chromosomes() { return mChromosomes; }

    public String chrPair() { return formChromosomePair(mChromosomes[SE_START], mChromosomes[SE_END]); }

    public static boolean validPositions(final long[] position) { return position[SE_START] > 0 && position[SE_END] > 0; }

}
