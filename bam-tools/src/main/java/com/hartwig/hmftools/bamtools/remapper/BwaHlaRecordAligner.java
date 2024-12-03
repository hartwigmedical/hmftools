package com.hartwig.hmftools.bamtools.remapper;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.assembly.alignment.Aligner;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class BwaHlaRecordAligner implements HlaRecordAligner
{

    private @NotNull
    final Aligner aligner;

    public BwaHlaRecordAligner(@NotNull final Aligner aligner)
    {
        this.aligner = aligner;
    }

    @NotNull
    @Override
    public List<SAMRecord> alignPair(@NotNull final RecordPair pair)
    {
        ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignments =
                aligner.alignSequences(pair.left.getReadBases(), pair.right.getReadBases());

        if (pair.left.getReadName().equals("A00624:8:HHKYHDSXX:2:1516:12156:4225")) {
            System.out.println("--------");
        }
        List<SAMRecord> result = new ArrayList<>();
        final List<SAMRecord> realignedRecordsLeft = produceRealignments(pair.left, alignments.getLeft());
        final List<SAMRecord> realignedRecordsRight = produceRealignments(pair.right, alignments.getRight());
        // The realigned records may still have mate references that point to hla contigs.
        // These need to be adjusted using the alignment info from the first record in the
        // realigned records for the complementary record.
        // TODO what if either of the returned sequences is empty?
        SAMRecord firstRealignedLeftRecord = realignedRecordsLeft.get(0);
        SAMRecord firstRealignedRightRecord = realignedRecordsRight.get(0);
        realignedRecordsLeft.forEach(record -> {
            record.setMateReferenceIndex(firstRealignedRightRecord.getReferenceIndex());
            record.setMateAlignmentStart(firstRealignedRightRecord.getAlignmentStart());
            result.add(record);
        });
        realignedRecordsRight.forEach(record -> {
            record.setMateReferenceIndex(firstRealignedLeftRecord.getReferenceIndex());
            record.setMateAlignmentStart(firstRealignedLeftRecord.getAlignmentStart());
            result.add(record);
        });
        return result;
    }

    @NotNull
    @Override
    public List<SAMRecord> alignRecord(@NotNull final SAMRecord record)
    {
        return produceRealignments(record, aligner.alignSequence(record.getReadBases()));
    }

    @NotNull
    private List<SAMRecord> produceRealignments(@NotNull final SAMRecord record, @NotNull final List<BwaMemAlignment> alignments)
    {
        return alignments
                .stream()
                .map(a -> HlaRecordAligner.createRemappedRecord(record, a))
                .collect(Collectors.toList());
    }
}
