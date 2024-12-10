package com.hartwig.hmftools.bamtools.remapper;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.esvee.assembly.alignment.Aligner;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class BwaHlaRecordAligner implements HlaRecordAligner
{

    private @NotNull
    final Aligner aligner;
    private final SAMFileHeader header;

    public BwaHlaRecordAligner(@NotNull final Aligner aligner, final SAMFileHeader newHeader)
    {
        this.aligner = aligner;
        header = newHeader;
    }

    @NotNull
    @Override
    public List<SAMRecord> alignPair(@NotNull final RecordPair pair)
    {
        ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignments =
                aligner.alignSequences(pair.leftBasesForRealignment(), pair.rightBasesForRealignment());

        List<SAMRecord> result = new ArrayList<>();
        final List<SAMRecord> realignedRecordsLeft = produceRealignments(pair.first, alignments.getLeft());
        final List<SAMRecord> realignedRecordsRight = produceRealignments(pair.second, alignments.getRight());
        // The realigned records may still have mate references that point to hla contigs.
        // These need to be adjusted using the alignment info from the first record in the
        // realigned records for the complementary record.
        // TODO what if either of the returned sequences is empty?
        SAMRecord firstRealignedLeftRecord = realignedRecordsLeft.stream().filter(sr -> !sr.isSecondaryOrSupplementary()).findFirst().get();
        SAMRecord firstRealignedRightRecord =
                realignedRecordsRight.stream().filter(sr -> !sr.isSecondaryOrSupplementary()).findFirst().get();
        if(pair.first.getReadName().equals("A00624:8:HHKYHDSXX:1:1664:13205:2472"))
        {
            System.out.println(pair);
        }
        if(pair.first.getReadName().equals("A00624:8:HHKYHDSXX:3:1536:28971:18630"))
        {
            System.out.println(pair);
        }
        fixProperPairFlag(firstRealignedLeftRecord, pair);
        fixProperPairFlag(firstRealignedRightRecord, pair);

        realignedRecordsLeft.forEach(record ->
        {
            setMateProperties(record, firstRealignedRightRecord);
            result.add(record);
        });
        realignedRecordsRight.forEach(record ->
        {
            setMateProperties(record, firstRealignedLeftRecord);
            result.add(record);
        });
        return result;
    }

    private static void fixProperPairFlag(SAMRecord record, RecordPair pair)
    {
        record.setProperPairFlag(pair.isProperPair());
    }

    private static void setMateProperties(SAMRecord record, SAMRecord mate)
    {
        record.setMateReferenceIndex(mate.getReferenceIndex());
        record.setMateAlignmentStart(mate.getAlignmentStart());
        record.setInferredInsertSize(calculateInsertSize(record, mate));
        record.setAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE, mate.getCigarString());
        record.setAttribute(SamRecordUtils.MATE_QUALITY_ATTRIBUTE, mate.getMappingQuality());
    }

    private static int calculateInsertSize(@NotNull final SAMRecord record, @NotNull final SAMRecord mate)
    {
        if(record.isSecondaryOrSupplementary())
        {
            return 0;
        }
        if(!Objects.equals(record.getReferenceIndex(), mate.getReferenceIndex()))
        {
            return 0;
        }

        if (record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
        {
            return 0;
        }
        if(record.getReadNegativeStrandFlag())
        {
            return -1 * (record.getAlignmentEnd() - record.getMateAlignmentStart() + 1);
        }
//        int cigarLength = CigarUtils.cigarBaseLength(CigarUtils.cigarFromStr(record.getCigarString()));
        return mate.getAlignmentEnd() - record.getAlignmentStart() + 1;
    }

    @NotNull
    private List<SAMRecord> produceRealignments(@NotNull final SAMRecord record, @NotNull final List<BwaMemAlignment> alignments)
    {
        return alignments
                .stream()
                .map(a -> remappedRecord(record, a))
                .collect(Collectors.toList());
    }

    private SAMRecord remappedRecord(@NotNull final SAMRecord original, @NotNull final BwaMemAlignment alignment)
    {
        RawFastaData raw = RawFastaData.fromRecord(original);
        SAMRecord remappedRecord = new SAMRecord(header);
        remappedRecord.setReadName(original.getReadName());
        remappedRecord.setHeader(this.header);
        remappedRecord.setFlags(alignment.getSamFlag());
        remappedRecord.setReferenceIndex(alignment.getRefId());
        remappedRecord.setAlignmentStart(alignment.getRefStart() + 1); // BwaMemAlignment is 0-based
        remappedRecord.setMappingQuality(alignment.getMapQual());
        remappedRecord.setCigarString(alignment.getCigar());
        if(remappedRecord.getReadNegativeStrandFlag())
        {
            remappedRecord.setReadBases(Nucleotides.reverseComplementBases(raw.bases));
            remappedRecord.setBaseQualities(Arrays.reverseArray(raw.qualities));
        }
        else
        {
            remappedRecord.setReadBases(raw.bases);
            remappedRecord.setBaseQualities(raw.qualities);
        }
        remappedRecord.setAttribute(SamRecordUtils.NUM_MUTATONS_ATTRIBUTE, alignment.getNMismatches());
        remappedRecord.setAttribute(SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE, alignment.getMDTag());
        remappedRecord.setAttribute(SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE, alignment.getAlignerScore());
        remappedRecord.setAttribute(SamRecordUtils.SUBOPTIMAL_SCORE_ATTRIBUTE, alignment.getSuboptimalScore());
        return remappedRecord;
    }
}
