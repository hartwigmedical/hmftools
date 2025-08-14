package com.hartwig.hmftools.redux.jitter;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.NONE;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.common.sequencing.BiomodalBamUtils.LOW_QUAL_CUTOFF;
import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.common.sequencing.SequencingType;

import org.apache.commons.lang3.Validate;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public abstract class ConsensusMarker
{
    public abstract ConsensusType consensusType(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record);

    @Nullable
    public static ConsensusMarker create(final JitterAnalyserConfig config)
    {
        SequencingType sequencingType = config.Sequencing;
        if(sequencingType == ILLUMINA && config.UsesDuplexUMIs)
            return new IlluminaDuplexUMIsConsensusMarker();

        if(sequencingType == SBX)
            return new SBXConsensusMarker();

        if(sequencingType == BIOMODAL)
            return new BiomodalConsensusMarker();

        return null;
    }

    private static final int INVALID_INDEX = -1;

    private static Pair<Integer, Integer> getMicrosatelliteBoundaries(
            final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record)
    {
        List<CigarElement> cigarElements = record.getCigar().getCigarElements();
        int cigarIdx = 0;
        int readIdx = 0;
        if(cigarElements.get(0).getOperator() == S)
        {
            cigarIdx++;
            readIdx += cigarElements.get(0).getLength();
        }
        else if(cigarElements.get(0).getOperator() == H)
        {
            cigarIdx++;
        }

        int refPos = record.getAlignmentStart();
        int startReadIdx = INVALID_INDEX;
        int endReadIdx = INVALID_INDEX;
        for(; cigarIdx < cigarElements.size(); cigarIdx++)
        {
            CigarElement el = cigarElements.get(cigarIdx);
            if(el.getOperator().isClipping())
                continue;

            boolean isRead = el.getOperator().consumesReadBases();
            boolean isRef = el.getOperator().consumesReferenceBases();

            if(isRead && isRef)
            {
                int endRefPos = refPos + el.getLength() - 1;
                if(startReadIdx == INVALID_INDEX && endRefPos >= refGenomeMicrosatellite.referenceStart())
                {
                    if(refPos >= refGenomeMicrosatellite.referenceStart())
                    {
                        startReadIdx = readIdx;
                    }
                    else
                    {
                        startReadIdx = readIdx + refGenomeMicrosatellite.referenceStart() - refPos;
                    }
                }

                if(endRefPos <= refGenomeMicrosatellite.referenceEnd())
                {
                    endReadIdx = readIdx + el.getLength() - 1;
                }
                else if(refPos <= refGenomeMicrosatellite.referenceEnd())
                {
                    endReadIdx = readIdx + refGenomeMicrosatellite.referenceEnd() - refPos;
                }

                readIdx += el.getLength();
                refPos += el.getLength();
            }
            else if(isRead)
            {
                if(startReadIdx == INVALID_INDEX && refPos - 1 >= refGenomeMicrosatellite.referenceStart())
                    startReadIdx = readIdx;

                if(refPos - 1 <= refGenomeMicrosatellite.referenceEnd())
                    endReadIdx = readIdx + el.getLength() - 1;

                readIdx += el.getLength();
            }
            else if(isRef)
            {
                refPos += el.getLength();
            }
            else
            {
                throw new IllegalStateException("Unreachable");
            }
        }

        Validate.isTrue(startReadIdx != INVALID_INDEX);
        Validate.isTrue(endReadIdx != INVALID_INDEX);

        return Pair.of(startReadIdx - 1, endReadIdx + 1);
    }

    public static class IlluminaDuplexUMIsConsensusMarker extends ConsensusMarker
    {
        @Override
        public ConsensusType consensusType(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record)
        {
            return extractConsensusType(record);
        }
    }

    public static class SBXConsensusMarker extends ConsensusMarker
    {
        @Override
        public ConsensusType consensusType(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record)
        {
            ConsensusType consensusType = extractConsensusType(record);

            if(consensusType != DUAL)
                return consensusType;

            int duplexBaseIndex = SbxBamUtils.extractDuplexBaseIndex(record);

            Pair<Integer,Integer> boundaries = getMicrosatelliteBoundaries(refGenomeMicrosatellite, record);

            if(record.getReadNegativeStrandFlag())
            {
                int endIndex = boundaries.getRight();
                return SbxBamUtils.inDuplexRegion(false, duplexBaseIndex, endIndex) ? DUAL : SINGLE;
            }
            else
            {
                int startIndex = boundaries.getLeft();
                return SbxBamUtils.inDuplexRegion(true, duplexBaseIndex, startIndex) ? DUAL : SINGLE;
            }
        }
    }

    public static class BiomodalConsensusMarker extends ConsensusMarker
    {
        @Override
        public ConsensusType consensusType(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record)
        {
            byte[] quals = record.getBaseQualities();
            Pair<Integer, Integer> boundaries = getMicrosatelliteBoundaries(refGenomeMicrosatellite, record);
            int startIdx = boundaries.getLeft();
            int endIdx = boundaries.getRight();
            int minQual = Integer.MAX_VALUE;
            for(int i = startIdx; i <= endIdx; i++)
                minQual = min(minQual, quals[i]);

            if(minQual > LOW_QUAL_CUTOFF)
                return DUAL;

            return NONE;
        }
    }
}
