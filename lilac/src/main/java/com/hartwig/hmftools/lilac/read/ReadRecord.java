package com.hartwig.hmftools.lilac.read;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.generateMappedCoords;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.List;
import java.util.stream.Collectors;

public class ReadRecord
{
    public final String Id;
    public final int SoftClippedStart; // soft-clipped bases at start and end
    public final int SoftClippedEnd;
    public final int PositionStart; // adjusted by soft-clipped bases
    public final int PositionEnd;
    public final int ReadStart;
    public final int ReadEnd;

    private final List<Indel> mIndels;
    private final SAMRecord mSamRecord;

    public ReadRecord(
            final String id, int softClippedStart, int softClippedEnd, final List<Indel> indels, int positionStart,
            int positionEnd, int readStart, int readEnd, final SAMRecord record)
    {
        Id = id;
        SoftClippedStart = softClippedStart;
        SoftClippedEnd = softClippedEnd;
        mIndels = indels;
        PositionStart = positionStart;
        PositionEnd = positionEnd;
        ReadStart = readStart;
        ReadEnd = readEnd;
        mSamRecord = record;
    }

    public String readInfo()
    {
        return format("%s:%d-%d %s",
                mSamRecord.getContig(), mSamRecord.getStart(), mSamRecord.getEnd(), mSamRecord.getCigarString());
    }

    public String toString() { return format("%s %s indels(%s) sc(%d/%d)",
            Id, readInfo(), mIndels.size(), SoftClippedStart, SoftClippedEnd); }

    public List<Indel> getIndels() { return mIndels; }

    public SAMRecord getSamRecord() { return mSamRecord; }

    public final int maxIndelSize()
    {
        return mIndels.stream().mapToInt(x -> abs(x.Length)).max().orElse(0);
    }

    public final boolean containsSoftClip() { return SoftClippedStart > 0 || SoftClippedEnd > 0; }

    public final boolean containsIndel() { return !mIndels.isEmpty(); }

    public final char[] codingRegionRead(boolean reverseCompliment)
    {
        final char[] readBases = forwardRead();

        if(!reverseCompliment)
            return readBases;

        final char[] reverseBases = new char[readBases.length];

        int index = 0;
        for(int i = readBases.length - 1; i >= 0; --i, ++index)
        {
            reverseBases[index] = reverseCompliment(readBases[i]);
        }

        return reverseBases;
    }

    public final int[] codingRegionQuality(boolean reverseCompliment)
    {
        final int[] readQuals = forwardQuality();

        if(!reverseCompliment)
            return readQuals;

        final int[] reverseQuals = new int[readQuals.length];

        int index = 0;
        for(int i = readQuals.length - 1; i >= 0; --i, ++index)
        {
            reverseQuals[index] = readQuals[i];
        }

        return reverseQuals;
    }

    private char[] forwardRead()
    {
        int readLength = ReadEnd - ReadStart + 1;
        final char[] readBases = new char[readLength];

        int index = 0;
        for(int i = ReadStart; i <= ReadEnd; ++i)
        {
            readBases[index++] = mSamRecord.getReadString().charAt(i);
        }

        return readBases;
    }

    private int[] forwardQuality()
    {
        int readLength = ReadEnd - ReadStart + 1;
        final int[] readQuals = new int[readLength];

        int index = 0;
        for(int i = ReadStart; i <= ReadEnd; ++i)
        {
            readQuals[index++] = mSamRecord.getBaseQualities()[i];
        }

        return readQuals;
    }

    public List<ReadRecord> alignmentsOnly()
    {
        final String chromosome = mSamRecord.getContig();
        final ChrBaseRegion outerRegion = new ChrBaseRegion(chromosome, PositionStart, PositionEnd);

        return mSamRecord.getAlignmentBlocks().stream()
                .map(x -> new ChrBaseRegion(chromosome, x.getReferenceStart(), x.getReferenceStart() + x.getLength()))
                .filter(x -> outerRegion.overlaps(x))
                .map(x -> new BaseRegion(max(outerRegion.start(), x.start()), min(outerRegion.end(), x.end())))
                .map(x -> create(x, mSamRecord, false, false))
                .collect(Collectors.toList());
    }

    public static ReadRecord create(final BaseRegion codingRegion, final SAMRecord record, boolean includeSoftClips, boolean includeIndels)
    {
        int softClipStart = leftSoftClipLength(record);
        int softClipEnd = rightSoftClipLength(record);
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();
        int recordStart = alignmentStart - softClipStart;
        int recordEnd = alignmentEnd + softClipEnd;

        int positionStart = max(codingRegion.start(), alignmentStart);
        int positionEnd = min(codingRegion.end(), alignmentEnd);

        if(record.getCigar().containsOperator(CigarOperator.N))
        {
            // handle splits which do not align with exon boundaries (which the coding region represents)
            List<int[]> mappedCoords = generateMappedCoords(record.getCigar(), record.getAlignmentStart());

            for(int[] mappedCoord : mappedCoords)
            {
                if(mappedCoord[SE_END] < codingRegion.start())
                    continue;

                if(mappedCoord[SE_START] > codingRegion.end())
                    break;

                positionStart = max(max(positionStart, mappedCoord[SE_START]), codingRegion.start());
                positionEnd = min(min(positionEnd, mappedCoord[SE_END]), codingRegion.end());
                break;
            }
        }

        int readIndexStart = record.getReadPositionAtReferencePosition(positionStart, true) - 1;
        int readIndexEnd = record.getReadPositionAtReferencePosition(positionEnd, true) - 1;
        positionStart = record.getReferencePositionAtReadPosition(readIndexStart + 1);
        positionEnd = record.getReferencePositionAtReadPosition(readIndexEnd + 1);

        // Add soft clip start
        if(positionStart == alignmentStart && softClipStart > 0 && includeSoftClips)
        {
            int earliestStart = max(codingRegion.start(), recordStart);
            readIndexStart = readIndexStart - positionStart + earliestStart;
            positionStart = earliestStart;
        }

        // Add soft clip end
        if(positionEnd == alignmentEnd && softClipEnd > 0 && includeSoftClips)
        {
            int latestEnd = min(codingRegion.end(), recordEnd);
            readIndexEnd = readIndexEnd + latestEnd - positionEnd;
            positionEnd = latestEnd;
        }

        int softClippedStart = max(alignmentStart - positionStart, 0);
        int softClippedEnd = max(0, positionEnd - alignmentEnd);

        List<Indel> indels = includeIndels ? indels(positionStart, positionEnd, record) : Lists.newArrayList();

        return new ReadRecord(
                record.getReadName(), softClippedStart, softClippedEnd, indels,
                positionStart, positionEnd, readIndexStart, readIndexEnd, record);
    }

    private static List<Indel> indels(int startPosition, int endPosition, final SAMRecord record)
    {
        List<Indel> indels = Lists.newArrayList();

        CigarHandler cigarHandler = new CigarHandler()
        {
            public void handleInsert(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                if(startPosition <= refPosition && refPosition <= endPosition)
                {
                    char baseChar = record.getReadString().charAt(readIndex);
                    indels.add(new Indel(record.getContig(), refPosition, String.valueOf(baseChar),
                            record.getReadString().substring(readIndex, readIndex + element.getLength() + 1)));
                }
            }

            public void handleDelete(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {

                if(startPosition <= refPosition && refPosition <= endPosition)
                {
                    char baseChar = record.getReadString().charAt(readIndex);
                    String delBases = String.valueOf(baseChar);

                    for(int i = 0; i < element.getLength(); ++i)
                        delBases += "N";

                    indels.add(new Indel(record.getContig(), refPosition, delBases, String.valueOf(baseChar)));
                }
            }
        };

        CigarHandler.traverseCigar(record, cigarHandler);
        return indels;
    }

    public static char reverseCompliment(char base)
    {
        switch(base)
        {
            case 'G': return 'C';
            case 'A': return 'T';
            case 'T': return 'A';
            case 'C': return 'G';
        }
        return base;
    }
}
