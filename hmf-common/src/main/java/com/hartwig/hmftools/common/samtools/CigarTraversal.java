package com.hartwig.hmftools.common.samtools;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public final class CigarTraversal
{
    public static void traverseCigar(final SAMRecord record, final CigarHandler handler)
    {
        final Cigar cigar = record.getCigar();

        int readIndex = 0;
        int refBase = record.getAlignmentStart();

        for(int i = 0; i < cigar.numCigarElements(); i++)
        {
            final CigarElement e = cigar.getCigarElement(i);
            switch(e.getOperator())
            {
                case H:
                    break; // ignore hard clips
                case P:
                    break; // ignore pads
                case S:
                    if(i == 0)
                    {
                        handler.handleLeftSoftClip(record, e);
                    }
                    else if(i == cigar.numCigarElements() - 1)
                    {
                        handler.handleRightSoftClip(record, e, readIndex, refBase);
                    }
                    readIndex += e.getLength();
                    break; // soft clip read bases
                case N:
                    handler.handleSkippedReference(record, e, readIndex - 1, refBase - 1);
                    refBase += e.getLength();
                    break;  // reference skip
                case D:
                    handler.handleDelete(record, e, readIndex - 1, refBase - 1);
                    refBase += e.getLength();
                    break;
                case I:
                    // TODO: Handle 1I150M
                    int refIndex = refBase - 1 - record.getAlignmentStart();
                    if(refIndex >= 0)
                    {
                        handler.handleInsert(record, e, readIndex - 1, refBase - 1);
                    }
                    readIndex += e.getLength();
                    break;
                case M:
                case EQ:
                case X:
                    boolean isFollowedByIndel = i < cigar.numCigarElements() - 1 && cigar.getCigarElement(i + 1).getOperator().isIndel();
                    final CigarElement element = isFollowedByIndel ? new CigarElement(e.getLength() - 1, e.getOperator()) : e;
                    handler.handleAlignment(record, element, readIndex, refBase);
                    readIndex += e.getLength();
                    refBase += e.getLength();
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with op: " + e.getOperator() + "in CIGAR: " + cigar);
            }
        }
    }

    public static int calcCigarLength(final String cigarStr)
    {
        int index = 0;
        int baseLength = 0;
        String basesStr = "";
        while(index < cigarStr.length())
        {
            char c = cigarStr.charAt(index);
            boolean isAddItem = (c == 'D' || c == 'M');
            boolean isIgnoreItem = (c == 'I' || c == 'N' || c == 'S' || c == 'H' || c == 'P' || c == '=' || c == 'X');

            if(isAddItem)
            {
                try { baseLength += Integer.parseInt(basesStr); } catch (Exception e) {}
                basesStr = "";
            }
            else if(isIgnoreItem)
            {
                basesStr = "";
            }
            else
            {
                basesStr += c;
            }

            ++index;
        }

        return baseLength;
    }
}
