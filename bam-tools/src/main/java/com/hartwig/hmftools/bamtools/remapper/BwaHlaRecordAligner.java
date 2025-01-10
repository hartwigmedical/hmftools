package com.hartwig.hmftools.bamtools.remapper;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
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

        // Get an alignments pair that minimises the distance between the read
        // and its mate.
        AlignmentsSelector alignmentsSelector = new AlignmentsSelector(leftAlignments, rightAlignments);
        HlaAlignmentPair bestAlignedPair = alignmentsSelector.closestAlignmentPair(refGenomeVersion);
        // Each element of the pair needs the other for its mate properties.
        SAMRecord principalLeftRemapped = bestAlignedPair.Left.createSamRecord(header, pair.leftData(), bestAlignedPair.Right);
//        SAMRecord principalLeftRemapped = remappedRecord(pair.leftData(), bestAlignedPair.Left);
//        SAMRecord principalRightRemapped = remappedRecord(pair.rightData(), bestAlignedPair.Right);
        SAMRecord principalRightRemapped = bestAlignedPair.Right.createSamRecord(header, pair.rightData(), bestAlignedPair.Left);
//        setMateProperties(principalLeftRemapped, principalRightRemapped);
//        setMateProperties(principalRightRemapped, principalLeftRemapped);
        // Because our alignments have been calculated by BWA one pair at a time,
        // they are missing the "proper pair" flag, which is set based on statistical
        // properties of large batches. Put this back if the pair are close together.
        fixProperPairFlag(principalLeftRemapped);
        fixProperPairFlag(principalRightRemapped);

        // Create a result list and add the pair.
        List<SAMRecord> result = new ArrayList<>();
        result.add(principalLeftRemapped);
        result.add(principalRightRemapped);

        // Calculate and add any supplementary alignments. These have their mate
        // properties set from the principal mate alignment.
        //                    setMateProperties(samRecord, principalRightRemapped);
        leftAlignments.supplementaryAlignments()
                .map(hla -> hla.createSamRecord(header, pair.leftData(), bestAlignedPair.Right))
                .forEach(result::add);
        //                    setMateProperties(samRecord, principalLeftRemapped);
        rightAlignments.supplementaryAlignments()
                .map(hla -> hla.createSamRecord(header, pair.rightData(), bestAlignedPair.Left))
                .forEach(result::add);

        return result;
    }

    private static void fixProperPairFlag(SAMRecord record)
    {
        int insertLength = Math.abs(record.getInferredInsertSize());
        final boolean insertLengthOK = insertLength < 1000 && insertLength > 50;
        final boolean sameStrand = Objects.equals(record.getReferenceIndex(), record.getMateReferenceIndex());
        if(insertLengthOK && sameStrand)
        {
            record.setProperPairFlag(true);
        }
    }

//    private static void setMateProperties(SAMRecord record, SAMRecord mate)
//    {
//        record.setMateReferenceIndex(mate.getReferenceIndex());
//        record.setMateAlignmentStart(mate.getAlignmentStart());
//        record.setInferredInsertSize(calculateInsertSize(record, mate));
//        record.setAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE, mate.getCigarString());
//        record.setAttribute(SamRecordUtils.MATE_QUALITY_ATTRIBUTE, mate.getMappingQuality());
//    }
//
//    private static int calculateInsertSize(@NotNull final SAMRecord record, @NotNull final SAMRecord mate)
//    {
//        if(record.isSecondaryOrSupplementary())
//        {
//            return 0;
//        }
//        if(!Objects.equals(record.getReferenceIndex(), mate.getReferenceIndex()))
//        {
//            return 0;
//        }
//
//        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
//        {
//            return 0;
//        }
//        if(record.getReadNegativeStrandFlag())
//        {
//            return -1 * (record.getAlignmentEnd() - record.getMateAlignmentStart() + 1);
//        }
//        return mate.getAlignmentEnd() - record.getAlignmentStart() + 1;
//    }
//
//    private SAMRecord remappedRecord(@NotNull final RawFastaData raw, @NotNull final HlaAlignment alignment)
//    {
//        SAMRecord remappedRecord = new SAMRecord(header);
//        remappedRecord.setReadName(raw.ReadName);
//        remappedRecord.setHeader(this.header);
//        remappedRecord.setFlags(alignment.getSamFlag());
//        remappedRecord.setReferenceIndex(alignment.getRefId());
//        remappedRecord.setAlignmentStart(alignment.getRefStart());
//        remappedRecord.setMappingQuality(alignment.getMapQual());
//        remappedRecord.setCigarString(alignment.getCigar());
//        if(remappedRecord.getReadNegativeStrandFlag())
//        {
//            remappedRecord.setReadBases(Nucleotides.reverseComplementBases(raw.Bases));
//            remappedRecord.setBaseQualities(Arrays.reverseArray(raw.Qualities));
//        }
//        else
//        {
//            remappedRecord.setReadBases(raw.Bases);
//            remappedRecord.setBaseQualities(raw.Qualities);
//        }
//        remappedRecord.setAttribute(SamRecordUtils.NUM_MUTATONS_ATTRIBUTE, alignment.getNMismatches());
//        remappedRecord.setAttribute(SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE, alignment.getMDTag());
//        remappedRecord.setAttribute(SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE, alignment.getAlignerScore());
//        remappedRecord.setAttribute(SamRecordUtils.SUBOPTIMAL_SCORE_ATTRIBUTE, alignment.getSuboptimalScore());
//        return remappedRecord;
//    }
}
