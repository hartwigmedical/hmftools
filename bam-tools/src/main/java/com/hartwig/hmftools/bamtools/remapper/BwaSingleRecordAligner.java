package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.Arrays;

import org.umccr.java.hellbender.utils.bwa.BwaMemAligner;
import org.umccr.java.hellbender.utils.bwa.BwaMemAlignment;
import org.umccr.java.hellbender.utils.bwa.BwaMemIndex;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class BwaSingleRecordAligner implements SingleRecordAligner
{
    private final BwaMemAligner mAligner;
    private final SAMFileHeader mHeader;

    public BwaSingleRecordAligner(BwaMemIndex index, SAMFileHeader header)
    {
        mAligner = new BwaMemAligner(index);
        mHeader = header;
        mAligner.dontInferPairEndStats();
        mAligner.setBandwidthOption(31);
    }

    @Override
    public List<SAMRecord> alignSequence(SAMRecord original)
    {
        List<BwaMemAlignment> bwaAlignments = mAligner.alignSeqs(List.of(original.getReadBases())).get(0);
        return bwaAlignments.stream()
                .map(bwaAlignment -> convertToSAM(bwaAlignment, original))
                .collect(Collectors.toList());
    }

    private SAMRecord convertToSAM(BwaMemAlignment bwaAlignment, SAMRecord original)
    {
        SAMRecord remappedRecord = original.deepCopy();
        remappedRecord.setHeader(mHeader);
        remappedRecord.setFlags(bwaAlignment.getSamFlag());

        remappedRecord.setReferenceIndex(bwaAlignment.getRefId());
        remappedRecord.setAlignmentStart(bwaAlignment.getRefStart() + 1);
        remappedRecord.setCigarString(bwaAlignment.getCigar());
        remappedRecord.setMappingQuality(bwaAlignment.getMapQual());

        if(remappedRecord.getReadNegativeStrandFlag())
        {
            remappedRecord.setReadBases(Nucleotides.reverseComplementBases(original.getReadBases()));
            remappedRecord.setBaseQualities(Arrays.reverseArray(original.getBaseQualities()));
        }
        else
        {
            remappedRecord.setReadBases(original.getReadBases());
            remappedRecord.setBaseQualities(original.getBaseQualities());
        }
        remappedRecord.setAttribute(SamRecordUtils.NUM_MUTATONS_ATTRIBUTE, bwaAlignment.getNMismatches());
        remappedRecord.setAttribute(SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE, bwaAlignment.getMDTag());
        remappedRecord.setAttribute(SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE, bwaAlignment.getAlignerScore());
        remappedRecord.setAttribute(SamRecordUtils.XS_ATTRIBUTE, bwaAlignment.getSuboptimalScore());

        return remappedRecord;
    }
}
