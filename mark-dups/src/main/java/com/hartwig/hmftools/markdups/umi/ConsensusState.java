package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.UNSET;
import static com.hartwig.hmftools.markdups.umi.UmiConfig.READ_ID_DELIM;

import java.util.List;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ConsensusState
{
    public byte[] Bases;
    public byte[] BaseQualities;
    public List<CigarElement> CigarElements;

    public int MinUnclippedPosStart;
    public int MaxUnclippedPosEnd;
    public int MinAlignedPosStart;
    public int MaxAlignedPosEnd;

    private ConsensusOutcome mOutcome;

    public ConsensusState()
    {
        Bases = null;
        BaseQualities = null;
        CigarElements = Lists.newArrayList();

        MinUnclippedPosStart = 0;
        MaxUnclippedPosEnd = 0;
        MinAlignedPosStart = 0;
        MaxAlignedPosEnd = 0;

        mOutcome = UNSET;
    }

    public ConsensusOutcome outcome() { return mOutcome; }
    public void setOutcome(final ConsensusOutcome outcome) { mOutcome = outcome; }

    void setBaseLength(int baseLength)
    {
        Bases = new byte[baseLength];
        BaseQualities = new byte[baseLength];
    }

    void setBoundaries(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();
        int readEnd = read.getAlignmentEnd();
        int unclippedStart = read.getCigar().isLeftClipped() ? readStart - read.getCigar().getFirstCigarElement().getLength() : readStart;
        int unclippedEnd = read.getCigar().isRightClipped() ? readEnd + read.getCigar().getLastCigarElement().getLength() : readEnd;

        if(MinUnclippedPosStart == 0)
        {
            MinUnclippedPosStart = unclippedStart;
            MaxUnclippedPosEnd = unclippedEnd;
            MinAlignedPosStart = readStart;
            MaxAlignedPosEnd = readEnd;
        }
        else
        {
            MinUnclippedPosStart = min(unclippedStart, MinUnclippedPosStart);
            MaxUnclippedPosEnd = max(unclippedEnd, MaxUnclippedPosEnd);
            MinAlignedPosStart = min(readStart, MinAlignedPosStart);
            MaxAlignedPosEnd = max(readEnd, MaxAlignedPosEnd);
        }
    }

    public SAMRecord createConsensusRead(final SAMRecord initialRead, final String groupIdentifier)
    {
        SAMRecord record = new SAMRecord(initialRead.getHeader());

        record.setReadName(formReadId(initialRead.getReadName(), groupIdentifier));
        record.setReadBases(Bases);
        record.setBaseQualities(BaseQualities);
        record.setReferenceName(initialRead.getReferenceName());

        record.setAlignmentStart(MinAlignedPosStart);
        record.setCigar(new Cigar(CigarElements));

        if(initialRead.getMateReferenceIndex() >= 0)
        {
            record.setMateReferenceName(initialRead.getMateReferenceName());
            record.setMateAlignmentStart(initialRead.getMateAlignmentStart());
            record.setMateReferenceIndex(initialRead.getMateReferenceIndex());
            record.setReadPairedFlag(true);
            record.setProperPairFlag(true);
        }
        else
        {
            record.setReadPairedFlag(false);
            record.setProperPairFlag(false);
        }

        record.setFlags(initialRead.getFlags());
        record.setDuplicateReadFlag(false); // being the new primary

        if(initialRead.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, initialRead.getAttribute(SUPPLEMENTARY_ATTRIBUTE));
        }

        record.setInferredInsertSize(initialRead.getInferredInsertSize());

        return record;
    }

    private static String formReadId(final String templateReadId, final String groupIdentifier)
    {
        int lastDelim = templateReadId.lastIndexOf(READ_ID_DELIM);
        return lastDelim > 0 ? templateReadId.substring(0, lastDelim) + READ_ID_DELIM + "CNS_" + groupIdentifier
                : templateReadId + READ_ID_DELIM + "CNS_" + groupIdentifier;
    }
}
