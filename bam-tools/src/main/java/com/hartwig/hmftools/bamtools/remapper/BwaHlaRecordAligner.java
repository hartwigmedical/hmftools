package com.hartwig.hmftools.bamtools.remapper;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
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
    private @NotNull
    final SAMFileHeader header;
    private @NotNull
    final RefGenomeVersion refGenomeVersion = RefGenomeVersion.V38; // TODO make parameter

    public BwaHlaRecordAligner(@NotNull final Aligner aligner, @NotNull final SAMFileHeader newHeader)
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
        AlignmentsList leftAlignments = new AlignmentsList(alignments.getLeft());
        AlignmentsList rightAlignments = new AlignmentsList(alignments.getRight());
        AlignmentsSelector alignmentsSelector = new AlignmentsSelector(leftAlignments, rightAlignments);

        final List<SAMRecord> realignedRecordsLeft = createNewRecords(alignmentsSelector.preferredLeftAlignments(), pair.leftData());
        final List<SAMRecord> realignedRecordsRight = createNewRecords(alignmentsSelector.preferredRightAlignments(), pair.rightData());
        // The realigned records may still have mate references that point to hla contigs.
        // These need to be adjusted using the alignment info from the first record in the
        // realigned records for the complementary record.
        SAMRecord firstRealignedLeftRecord = realignedRecordsLeft.stream().filter(sr -> !sr.isSecondaryOrSupplementary()).findFirst().get();
        SAMRecord firstRealignedRightRecord =
                realignedRecordsRight.stream().filter(sr -> !sr.isSecondaryOrSupplementary()).findFirst().get();

        List<SAMRecord> result = new ArrayList<>();
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
        fixProperPairFlag(firstRealignedLeftRecord, pair);
        fixProperPairFlag(firstRealignedRightRecord, pair);
        return result;
    }

    private List<SAMRecord> createNewRecords(List<PreferredAlignment> alignments, RawFastaData data)
    {
        return alignments.stream().map(preferredAlignment -> remappedRecord(data, preferredAlignment))
                .collect(Collectors.toList());
    }

    private static void fixProperPairFlag(SAMRecord record, RecordPair pair)
    {
        int insertLength = Math.abs(record.getInferredInsertSize());
        if(insertLength < 1200 && insertLength > 50)
        {
            record.setProperPairFlag(pair.isProperPair());
            //        if(pair.isProperPair())
            //        {
            int newQuality = Math.min(60, record.getMappingQuality() + 20);
            record.setMappingQuality(newQuality);
            //        }
        }
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

        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
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
    private List<SAMRecord> produceRealignments(
            @NotNull final RawFastaData raw,
            @NotNull final List<BwaMemAlignment> alignments,
            @NotNull final BwaMemAlignment mate)
    {
        return alignments
                .stream()
                .map(a -> new PreferredAlignment(a, mate, refGenomeVersion))
                .map(a -> remappedRecord(raw, a))
                .collect(Collectors.toList());
    }

    private SAMRecord remappedRecord(@NotNull final RawFastaData raw, @NotNull final PreferredAlignment alignment)
    {
        SAMRecord remappedRecord = new SAMRecord(header);
        remappedRecord.setReadName(raw.readName);
        remappedRecord.setHeader(this.header);
        remappedRecord.setFlags(alignment.getSamFlag());
        remappedRecord.setReferenceIndex(alignment.getRefId());
        remappedRecord.setAlignmentStart(alignment.getRefStart());
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
