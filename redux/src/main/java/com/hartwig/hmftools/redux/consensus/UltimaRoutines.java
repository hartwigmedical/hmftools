package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.XS_ATTRIBUTE;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.HALF_PHRED_SCORE_SCALING;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_T0_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_TP_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG_DELIM;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG_INDEX_DELIM;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractLowQualCount;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractT0Values;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.isHighBaseQual;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

public final class UltimaRoutines
{
    private static final Set<String> NON_ULTIMA_REQUIRED_ATTRIBUTES = Set.of(
            SAMTag.AS.name(), SAMTag.RG.name(), SAMTag.NM.name(), XS_ATTRIBUTE, SUPPLEMENTARY_ATTRIBUTE);

    private static final Set<String> ULTIMA_RAW_QUAL_ATTRIBUTES = Set.of(ULTIMA_T0_TAG, ULTIMA_TP_TAG);

    private static boolean UTLIMA_DROP_QUAL_TAGS = false;
    private static final String CFG_ULTIMA_DROP_QUA_TAGS = "ultima_drop_qual_tags";

    public static void registerUltimaConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(CFG_ULTIMA_DROP_QUA_TAGS, "Drop Ultima qual tags");
    }

    public static void setUltimaConfig(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasFlag(CFG_ULTIMA_DROP_QUA_TAGS))
            UTLIMA_DROP_QUAL_TAGS = true;
    }

    public static void stripAttributes(final SAMRecord record, @Nullable final UltimaStats stats)
    {
        boolean keepUltimaRequired = !record.getReadUnmappedFlag() && !record.isSecondaryAlignment()
                && !record.getSupplementaryAlignmentFlag() && !record.getDuplicateReadFlag();

        if(UTLIMA_DROP_QUAL_TAGS)
            keepUltimaRequired = false;

        boolean formCustomTags = !record.getReadUnmappedFlag() && !record.isSecondaryAlignment();

        stripAttributes(record, formCustomTags, keepUltimaRequired);

        if(formCustomTags && stats != null)
        {
            int lowQualBaseCount = extractLowQualCount(record);
            stats.addLowQualCount(lowQualBaseCount);
        }
    }

    public static void stripAttributes(final SAMRecord record, boolean formCustomTags, boolean keepUltimaRawQuals)
    {
        ConsensusType consensusType = ConsensusType.NONE;
        String lowQualTag = null;

        if(formCustomTags)
        {
            consensusType = UltimaBamUtils.deriveConsensusType(record);
            lowQualTag = formLowQualTag(record);
        }

        // clean-up unrequired attributes
        List<SAMRecord.SAMTagAndValue> existingAttributes = record.getAttributes();

        record.clearAttributes();

        for(SAMRecord.SAMTagAndValue attribute : existingAttributes)
        {
            if(NON_ULTIMA_REQUIRED_ATTRIBUTES.contains(attribute.tag))
                record.setAttribute(attribute.tag, attribute.value);

            if(keepUltimaRawQuals && ULTIMA_RAW_QUAL_ATTRIBUTES.contains(attribute.tag))
                record.setAttribute(attribute.tag, attribute.value);
        }

        if(formCustomTags)
        {
            record.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, consensusType.toString());
            record.setAttribute(ULT_QUAL_TAG, lowQualTag);
        }
    }
    public static void preProcessRead(final SAMRecord record, final UltimaStats stats)
    {
        stripAttributes(record, stats);
    }

    public static void finaliseRead(final RefGenomeInterface refGenome, final SAMRecord record)
    {
    /*
        ConsensusType consensusType = UltimaBamUtils.deriveConsensusType(record);
        record.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, consensusType.toString());
        setLowQualTag(record);
    */
    }

    public static String formLowQualTag(final SAMRecord record)
    {
        // tag has the form: 5=2,4-6,8 where
        // 5 is the total number of low-qual bases, followed by the specific indices for these - single index values and index ranges
        // these indices are determined from from homolpolymer complete deletions (t0 tag) and adjustments (tP tag)
        byte[] t0Values = extractT0Values(record);

        byte[] bases = record.getReadBases();
        byte[] baseQuals = record.getBaseQualities();

        StringJoiner sj = null;

        int hpStartIndex = -1;
        int hpEndIndex = -1;
        boolean hpLowQual = false;
        int lowQualStartIndex = -1;
        int totalLowQualCount = 0;

        for(int i = 0; i < baseQuals.length; ++i)
        {
            byte qual = baseQuals[i];

            boolean isLowQual = false;

            // deletions
            if(t0Values[i] < ULTIMA_HP_DELETION_LOW_QUAL_THRESHOLD)
                isLowQual = true;

            // homopolymer adjustments
            if(hpStartIndex >= 0 && i <= hpEndIndex)
            {
                // in an existing homopolymer with low qual
                isLowQual |= hpLowQual;
            }
            else
            {
                // assess a new homopolymer
                hpLowQual = false;

                if(!isHighBaseQual(qual))
                {
                    byte base = bases[i];

                    // find the length of the HP
                    int hpLength = 1;

                    for(int j = i + 1; j < baseQuals.length; ++j)
                    {
                        if(bases[j] != base)
                            break;

                        ++hpLength;
                    }

                    byte tpQual = qual;

                    if(hpLength > 1)
                    {
                        if(tpQual - HALF_PHRED_SCORE_SCALING < 0)
                            tpQual = 0;
                        else
                            tpQual -= HALF_PHRED_SCORE_SCALING;
                    }

                    hpStartIndex = i;
                    hpEndIndex = i + hpLength - 1;

                    if(belowLowQualThreshold(hpLength, tpQual))
                    {
                        isLowQual = true;
                        hpLowQual = true;
                    }
                }
            }

            if(isLowQual)
            {
                if(lowQualStartIndex < 0)
                    lowQualStartIndex = i;

                if(i == baseQuals.length - 1)
                {
                    if(sj == null)
                        sj = new StringJoiner(ULT_QUAL_TAG_INDEX_DELIM);

                    addLowQualTagEntry(sj, lowQualStartIndex, baseQuals.length - 1);
                    totalLowQualCount += baseQuals.length - lowQualStartIndex;
                }
            }
            else
            {
                if(lowQualStartIndex >= 0)
                {
                    int lowQualEndIndex = i - 1;

                    if(sj == null)
                        sj = new StringJoiner(ULT_QUAL_TAG_INDEX_DELIM);

                    addLowQualTagEntry(sj, lowQualStartIndex, lowQualEndIndex);
                    totalLowQualCount += lowQualEndIndex - lowQualStartIndex + 1;

                    lowQualStartIndex = -1; // reset
                }
            }
        }

        return sj != null ? format("%d=%s", totalLowQualCount, sj) : null;
    }

    private static void addLowQualTagEntry(final StringJoiner sj, int lowQualStartIndex, int lowQualEndIndex)
    {
        if(lowQualEndIndex > lowQualStartIndex)
            sj.add(format("%d-%d", lowQualStartIndex, lowQualEndIndex));
        else
            sj.add(String.valueOf(lowQualStartIndex));
    }

    private static final byte ULTIMA_HP_DELETION_LOW_QUAL_THRESHOLD = 20;

    private static boolean belowLowQualThreshold(int homopolymerLength, byte qual)
    {
        if(homopolymerLength == 1)
            return qual < 25;
        else if(homopolymerLength == 2)
            return qual < 22;
        else if(homopolymerLength == 3)
            return qual < 19;
        else if(homopolymerLength <= 6)
            return qual < 16;
        else if(homopolymerLength < 8)
            return qual < 14;
        else
            return qual < 12;
    }
}
