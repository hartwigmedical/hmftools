package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CompareUtils
{
    private static final List<String> KEY_ATTRIBUTES = List.of(SUPPLEMENTARY_ATTRIBUTE, MATE_CIGAR_ATTRIBUTE);

    @VisibleForTesting
    public static List<String> compareReads(final SAMRecord origRead, final SAMRecord newRead, final CompareConfig config)
    {
        List<String> diffs = new ArrayList<>();

        boolean hasUnmappedDifference = origRead.getReadUnmappedFlag() != newRead.getReadUnmappedFlag();
        boolean hasMateUnmappedDifference = origRead.getMateUnmappedFlag() != newRead.getMateUnmappedFlag();

        boolean skipUnmappingDifference = config.IgnoreReduxUnmapped && (hasUnmappedDifference || hasMateUnmappedDifference);

        if(!skipUnmappingDifference)
        {
            if(origRead.getInferredInsertSize() != newRead.getInferredInsertSize())
                diffs.add(format("insertSize(%d/%d)", origRead.getInferredInsertSize(), newRead.getInferredInsertSize()));

            if(!origRead.getCigarString().equals(newRead.getCigarString()))
                diffs.add(format("cigar(%s/%s)", origRead.getCigarString(), newRead.getCigarString()));

            if(origRead.getMappingQuality() != newRead.getMappingQuality())
                diffs.add(format("mapQuality(%d/%d)", origRead.getMappingQuality(), newRead.getMappingQuality()));
        }

        if(origRead.getFlags() != newRead.getFlags())
        {
            if(origRead.getReadNegativeStrandFlag() != newRead.getReadNegativeStrandFlag())
                diffs.add(format("negStrand(%s/%s)", origRead.getReadNegativeStrandFlag(), newRead.getReadNegativeStrandFlag()));

            if(!config.IgnoreDupDiffs && origRead.getDuplicateReadFlag() != newRead.getDuplicateReadFlag())
                diffs.add(format("duplicate(%s/%s)", origRead.getDuplicateReadFlag(), newRead.getDuplicateReadFlag()));

            if(!skipUnmappingDifference)
            {
                if(hasUnmappedDifference)
                    diffs.add(format("unmapped(%s/%s)", origRead.getReadUnmappedFlag(), newRead.getReadUnmappedFlag()));

                if(hasMateUnmappedDifference)
                    diffs.add(format("mateUnmapped(%s/%s)", origRead.getMateUnmappedFlag(), newRead.getMateUnmappedFlag()));
            }
        }

        // check key attributes:
        for(String attribute : KEY_ATTRIBUTES)
        {
            if(attribute.equals(SUPPLEMENTARY_ATTRIBUTE))
            {
                if(config.IgnoreSupplementaryAttribute)
                    continue;

                if(origRead.hasAttribute(UNMAP_ATTRIBUTE) || newRead.hasAttribute(UNMAP_ATTRIBUTE))
                    continue;
            }

            if(skipUnmappingDifference && attribute.equals(MATE_CIGAR_ATTRIBUTE))
                continue;

            String readAttr1 = origRead.getStringAttribute(attribute);
            String readAttr2 = newRead.getStringAttribute(attribute);

            if(!Objects.equals(readAttr1, readAttr2)) // Objects.equals handles null case
            {
                diffs.add(format("attrib_%s(%s/%s)", attribute,
                        readAttr1 == null ? "missing" : readAttr1,
                        readAttr2 == null ? "missing" : readAttr2));
            }
        }

        if(config.CheckBasesAndQuals)
        {
            if(!origRead.getReadString().equals(newRead.getReadString()))
            {
                diffs.add(format("bases(%s/%s)", origRead.getReadString(), newRead.getReadString()));
            }

            if(!origRead.getBaseQualityString().equals(newRead.getBaseQualityString()))
            {
                diffs.add(format("baseQual(%s/%s)", origRead.getBaseQualityString(), newRead.getBaseQualityString()));
            }
        }

        return diffs;
    }

    public static SamReaderFactory makeSamReaderFactory(CompareConfig config)
    {
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT);
        if(config.RefGenomeFile != null && !config.RefGenomeFile.isEmpty())
        {
            readerFactory.referenceSequence(new File(config.RefGenomeFile));
        }
        return readerFactory;
    }
}
