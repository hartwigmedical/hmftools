package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_MATCH;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.SPLICED_BOTH;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.SPLICED_ONE;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.UNKNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.lowerChromosome;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

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
                sjPosition = read.getCoordsBoundary(true);
                sjOrientation = -1;
            }
            else
            {
                sjPosition = read.getCoordsBoundary(false);
                sjOrientation = 1;
            }

            if(sjCount == 0)
            {
                ++sjCount;
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

                break;
            }
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

    public final List<ReadRecord> getReads() { return mReads; }
    public FusionFragmentType type() { return mType; }
    public final String[] chromosomes() { return mChromosomes; }

    public final long[] splicePositions() { return mSjPositions; }
    public final byte[] spliceOrientations() { return mSjOrientations; }
    public boolean hasValidSpliceData() { return mSjValid[SE_START] && mSjValid[SE_END]; }

    public String chrPair() { return formChromosomePair(mChromosomes[SE_START], mChromosomes[SE_END]); }

    public static boolean validPositions(final long[] position) { return position[SE_START] > 0 && position[SE_END] > 0; }

    public void populateGeneCandidates(final List<List<String>> spliceGeneIds)
    {
        // each fragment supporting the splice junction will have the same set of candidate genes
        for(int se = SE_START; se <= SE_END; ++se)
        {
            for(final ReadRecord read : mReads)
            {
                if(!read.Chromosome.equals(mChromosomes[se]))
                    continue;

                for(Map.Entry<RegionMatchType,List<TransExonRef>> entry : read.getTransExonRefs().entrySet())
                {
                    if(entry.getKey() != EXON_BOUNDARY && entry.getKey() != EXON_MATCH)
                        continue;

                    if(read.getCoordsBoundary(true) == mSjPositions[se] || read.getCoordsBoundary(false) == mSjPositions[se])
                    {
                        for(TransExonRef transData : entry.getValue())
                        {
                            if(!spliceGeneIds.get(se).contains(transData.GeneId))
                                spliceGeneIds.get(se).add(transData.GeneId);
                        }
                    }
                }
            }
        }
    }
}
