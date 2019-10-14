package com.hartwig.hmftools.sage.cigar;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class CigarTraversal {

    public static void traverseCigar(@NotNull final SAMRecord record, @NotNull final CigarHandler handler) {
        final Cigar cigar = record.getCigar();

        int readIndex = 0;
        int refBase = record.getAlignmentStart();

        for (final CigarElement e : cigar.getCigarElements()) {
            switch (e.getOperator()) {
                case H:
                    break; // ignore hard clips
                case P:
                    break; // ignore pads
                case S:
                    readIndex += e.getLength();
                    break; // soft clip read bases
                case N:
                    refBase += e.getLength();
                    break;  // reference skip
                case D:
                    handler.handleDelete(record, e, readIndex - 1, refBase - 1);
                    refBase += e.getLength();
                    break;
                case I:
                    // TODO: Handle 1I150M
                    int refIndex = refBase - 1 - record.getAlignmentStart();
                    if (refIndex >= 0) {
                        handler.handleInsert(record, e, readIndex - 1, refBase - 1);
                    }
                    readIndex += e.getLength();
                    break;
                case M:
                case EQ:
                case X:
                    handler.handleAlignment(record, e, readIndex, refBase);
                    readIndex += e.getLength();
                    refBase += e.getLength();
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with  op: " + e.getOperator() + "in CIGAR: " + cigar);
            }
        }

    }

}
