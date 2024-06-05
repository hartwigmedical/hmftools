package com.hartwig.hmftools.common.bam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public interface CigarHandler
{
    default void handleLeftSoftClip(final SAMRecord record, final CigarElement element) {}

    default void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    default void handleAlignment(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    // Note that for insert, readIndex and refPosition are BEFORE the event
    default void handleInsert(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    // Note that for delete, readIndex and refPosition are BEFORE the event
    default void handleDelete(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    // Note that for skipped, readIndex and refPosition are BEFORE the event
    default void handleSkippedReference(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    static void traverseCigar(final SAMRecord record, final CigarHandler handler)
    {
        final Cigar cigar = record.getCigar();

        int readIndex = 0;
        int refBase = record.getAlignmentStart();

        for(int i = 0; i < cigar.numCigarElements(); i++)
        {
            final CigarElement element = cigar.getCigarElement(i);

            switch(element.getOperator())
            {
                case H: // ignore hard clips - no need to skip either bases or positions
                case P: // ignore pads
                    break;

                case S:
                    if(i == 0)
                    {
                        handler.handleLeftSoftClip(record, element);
                    }
                    else if(i == cigar.numCigarElements() - 1)
                    {
                        handler.handleRightSoftClip(record, element, readIndex, refBase);
                    }
                    readIndex += element.getLength();
                    break;

                case N:
                    handler.handleSkippedReference(record, element, readIndex - 1, refBase - 1);
                    refBase += element.getLength();
                    break;

                case D:
                    handler.handleDelete(record, element, readIndex - 1, refBase - 1);
                    refBase += element.getLength();
                    break;

                case I:

                    if(readIndex == 0) // unusual and ref base is still assumed to be the base prior
                        handler.handleInsert(record, element, 0, refBase - 1);
                    else
                        handler.handleInsert(record, element, readIndex - 1, refBase - 1);

                    readIndex += element.getLength();
                    break;

                case M:
                case EQ:
                case X:
                    handler.handleAlignment(record, element, readIndex, refBase);
                    readIndex += element.getLength();
                    refBase += element.getLength();
                    break;

                default:
                    throw new IllegalStateException("Case statement didn't deal with op: " + element.getOperator() + "in CIGAR: " + cigar);
            }
        }
    }

}
