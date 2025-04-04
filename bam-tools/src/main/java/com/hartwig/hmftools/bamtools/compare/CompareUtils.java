package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.codon.Nucleotides.complement;

import static htsjdk.samtools.util.SequenceUtil.reverseComplement;

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

        // check the read bases, make sure we account for the read negative strand flag
        if(!basesMatch(
                origRead.getReadString(), origRead.getReadNegativeStrandFlag(),
                newRead.getReadString(), newRead.getReadNegativeStrandFlag()))
        {
            diffs.add(format("bases(%s/%s)",
                    origRead.getReadNegativeStrandFlag() ? reverseComplement(origRead.getReadString()) : origRead.getReadString(),
                    newRead.getReadNegativeStrandFlag() ? reverseComplement(newRead.getReadString()) : newRead.getReadString()));
        }

        // check the base qual, make sure we account for the read negative strand flag
        if(!stringsMatch(
                origRead.getBaseQualityString(), origRead.getReadNegativeStrandFlag(),
                newRead.getBaseQualityString(), newRead.getReadNegativeStrandFlag()))
        {
            diffs.add(format("baseQual(%s/%s)",
                    origRead.getReadNegativeStrandFlag()
                            ? new StringBuilder(origRead.getBaseQualityString()).reverse()
                            : origRead.getBaseQualityString(),
                    newRead.getReadNegativeStrandFlag()
                            ? new StringBuilder(newRead.getBaseQualityString()).reverse()
                            : newRead.getBaseQualityString()));
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

    // check if the string match, with the extra caveat that they could be reversed
    public static boolean stringsMatch(String str1, boolean str1Reversed, String str2, boolean str2Reversed)
    {
        if(str1Reversed == str2Reversed)
        {
            return str1.equals(str2);
        }

        // they are not the same orientation, check str1 forward and str2 backwards
        if (str1.length() == str2.length())
        {
            for(int i = 0; i < str1.length(); i++)
            {
                if(str1.charAt(i) != str2.charAt(str2.length() - i - 1))
                {
                    return false;
                }
            }
        }
        return true;
    }

    // check if the bases match, with the extra caveat that they could be reverse complemented
    public static boolean basesMatch(String bases1, boolean bases1Reversed, String bases2, boolean bases2Reversed)
    {
        if(bases1Reversed == bases2Reversed)
        {
            return bases1.equals(bases2);
        }

        // they are not the same orientation, check bases1 forward and bases2 backwards and apply complement
        if(bases1.length() == bases2.length())
        {
            for(int i = 0; i < bases1.length(); i++)
            {
                // need to be complemented
                if(complement(bases1.charAt(i)) != bases2.charAt(bases2.length() - i - 1))
                {
                    return false;
                }
            }
        }
        return true;
    }
}
