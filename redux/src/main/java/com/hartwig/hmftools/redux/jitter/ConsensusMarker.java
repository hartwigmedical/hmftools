package com.hartwig.hmftools.redux.jitter;

import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractConsensusType;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.redux.ReduxConfig.SEQUENCING_TYPE;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public abstract class ConsensusMarker
{
    public abstract ConsensusType consensusType(final MicrosatelliteSite microsatelliteSite, final SAMRecord record);

    public static ConsensusMarker create()
    {
        if(SEQUENCING_TYPE == SBX)
            return new SBXConsensusMarker();

        return new StandardConsensusMarker();
    }

    private static final int INVALID_INDEX = -1;

    private static Pair<Integer, Integer> getMicrosatelliteBoundaries(final MicrosatelliteSite microsatelliteSite, final SAMRecord record)
    {
        List<CigarElement> cigarElements = record.getCigar().getCigarElements();
        int cigarIndex = 0;
        int readIndex = 0;
        if(cigarElements.get(0).getOperator() == S)
        {
            cigarIndex++;
            readIndex += cigarElements.get(0).getLength();
        }
        else if(cigarElements.get(0).getOperator() == H)
        {
            cigarIndex++;
        }

        int refPosition = record.getAlignmentStart();
        int startReadIndex = INVALID_INDEX;
        int endReadIndex = INVALID_INDEX;

        for(; cigarIndex < cigarElements.size(); cigarIndex++)
        {
            CigarElement el = cigarElements.get(cigarIndex);
            if(el.getOperator().isClipping())
                continue;

            boolean isRead = el.getOperator().consumesReadBases();
            boolean isRef = el.getOperator().consumesReferenceBases();

            if(isRead && isRef)
            {
                int endRefPos = refPosition + el.getLength() - 1;
                if(startReadIndex == INVALID_INDEX && endRefPos >= microsatelliteSite.referenceStart())
                {
                    if(refPosition >= microsatelliteSite.referenceStart())
                    {
                        startReadIndex = readIndex;
                    }
                    else
                    {
                        startReadIndex = readIndex + microsatelliteSite.referenceStart() - refPosition;
                    }
                }

                if(endRefPos <= microsatelliteSite.referenceEnd())
                {
                    endReadIndex = readIndex + el.getLength() - 1;
                }
                else if(refPosition <= microsatelliteSite.referenceEnd())
                {
                    endReadIndex = readIndex + microsatelliteSite.referenceEnd() - refPosition;
                }

                readIndex += el.getLength();
                refPosition += el.getLength();
            }
            else if(isRead)
            {
                if(startReadIndex == INVALID_INDEX && refPosition - 1 >= microsatelliteSite.referenceStart())
                    startReadIndex = readIndex;

                if(refPosition - 1 <= microsatelliteSite.referenceEnd())
                    endReadIndex = readIndex + el.getLength() - 1;

                readIndex += el.getLength();
            }
            else if(isRef)
            {
                refPosition += el.getLength();
            }
            else
            {
                throw new IllegalStateException("Unreachable");
            }
        }

        return Pair.of(startReadIndex - 1, endReadIndex + 1);
    }

    public static class StandardConsensusMarker extends ConsensusMarker
    {
        @Override
        public ConsensusType consensusType(final MicrosatelliteSite microsatelliteSite, final SAMRecord record)
        {
            return extractConsensusType(record);
        }
    }

    public static class SBXConsensusMarker extends ConsensusMarker
    {
        @Override
        public ConsensusType consensusType(final MicrosatelliteSite microsatelliteSite, final SAMRecord record)
        {
            ConsensusType consensusType = extractConsensusType(record);

            if(consensusType != DUAL)
                return consensusType;

            int duplexBaseIndex = SbxBamUtils.extractDuplexBaseIndex(record);

            Pair<Integer,Integer> boundaries = getMicrosatelliteBoundaries(microsatelliteSite, record);

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
}
