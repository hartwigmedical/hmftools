package com.hartwig.hmftools.svtools.sequence;

public class InsertSeqData
{
    public final String SampleId;
    public final int SvId;
    public final String VcfId;
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public final String InsertSeq;
    public final double CopyNumber;
    public final double CopyNumberChange;
    public final String LinkedBy;
    public final String RefContext;
    public final String Alignments;

    public InsertSeqData(
            final String sampleId, final int svId, final String vcfId, final String chromosome, final int position,
            final byte orientation, final String insertSeq, final double copyNumber, final double copyNumberChange, final String linkedBy,
            final String refContext, final String alignments)
    {
        SampleId = sampleId;
        SvId = svId;
        VcfId = vcfId;
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        InsertSeq = insertSeq;
        CopyNumber = copyNumber;
        CopyNumberChange = copyNumberChange;
        LinkedBy = linkedBy;
        RefContext = refContext;
        Alignments = alignments;
    }

    public static InsertSeqData fromCsv(final String line)
    {
        String[] values = line.split(",", -1);
        int index = 0;

        // SampleId,SvId,VcfId,StartChromosome,StartPosition,StartOrientation,InsertSequence,AdjustedCopyNumberStart,AdjustedCopyNumberChangeStart,
        // StartLinkedBy,StartRefContext,InsertSequenceAlignments,InsertSequenceRepeatType
        return new InsertSeqData(
                values[index++], Integer.parseInt(values[index++]), values[index++], values[index++], Integer.parseInt(values[index++]),
                Byte.parseByte(values[index++]), values[index++], Double.parseDouble(values[index++]), Double.parseDouble(values[index++]),
                values[index++], values[index++], values[index++]);
    }
}
