package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.Set;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;

public class RefGenomeCompatibility
{
    private static final String SAME_AS_REFERENCE_NAME = "=";

    private final Set<String> mAllowedContigs;
    private final SAMFileHeader mOutputHeader;

    public RefGenomeCompatibility(final Set<String> allowedContigs, final SAMFileHeader inputHeader, final SAMSequenceDictionary refDictionary)
    {
        mAllowedContigs = allowedContigs;
        mOutputHeader = inputHeader.clone();
        mOutputHeader.setSequenceDictionary(refDictionary);
    }

    public SAMFileHeader outputHeader()
    {
        return mOutputHeader.clone();
    }

    public boolean recordIsCompatible(final SAMRecord record)
    {
        if(!contigIsCompatible(record.getReferenceName()))
            return false;

        if(!contigIsCompatible(resolvedMateReferenceName(record)))
            return false;

        String suppData = record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE);
        if(suppData == null)
            return true;

        return SupplementaryReadData.extractAlignments(suppData).stream().allMatch(x -> contigIsCompatible(x.Chromosome));
    }

    public boolean prepareForWrite(final SAMRecord record)
    {
        if(!recordIsCompatible(record))
            return false;

        String referenceName = record.getReferenceName();
        String mateReferenceName = resolvedMateReferenceName(record);

        record.setHeader(mOutputHeader);
        record.setReferenceName(referenceName);
        record.setMateReferenceName(mateReferenceName);

        return true;
    }

    private boolean contigIsCompatible(final String contig)
    {
        return contig == null || contig.equals(NO_CHROMOSOME_NAME) || mAllowedContigs.contains(contig);
    }

    private static String resolvedMateReferenceName(final SAMRecord record)
    {
        String mateReferenceName = record.getMateReferenceName();
        return SAME_AS_REFERENCE_NAME.equals(mateReferenceName) ? record.getReferenceName() : mateReferenceName;
    }
}
