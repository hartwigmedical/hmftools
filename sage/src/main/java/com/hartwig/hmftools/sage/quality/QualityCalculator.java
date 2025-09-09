package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.minQual;
import static com.hartwig.hmftools.sage.ReferenceData.isHighlyPolymorphic;
import static com.hartwig.hmftools.sage.SageConfig.SEQUENCING_TYPE;
import static com.hartwig.hmftools.sage.SageConfig.isSbx;
import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.SageConstants.HIGHLY_POLYMORPHIC_GENES_MAX_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.MAX_MAP_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.READ_EDGE_PENALTY_0;
import static com.hartwig.hmftools.sage.SageConstants.READ_EDGE_PENALTY_1;
import static com.hartwig.hmftools.sage.seqtech.SbxUtils.indelCoreQuality;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.seqtech.SbxUtils;
import com.hartwig.hmftools.sage.seqtech.UltimaUtils;

import htsjdk.samtools.SAMRecord;

public class QualityCalculator
{
    private final QualityConfig mConfig;
    private final BqrRecordMap mQualityRecalibrationMap;
    private final MsiJitterCalcs mMsiJitterCalcs;
    private final RefSequence mRefBases;

    public static final byte INVALID_BASE_QUAL = -1;

    public QualityCalculator(
            final SageConfig config, final BqrRecordMap qualityRecalibrationMap, final RefSequence refBases,
            final MsiJitterCalcs msiJitterCalcs)
    {
        mConfig = config.Quality;
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mMsiJitterCalcs = msiJitterCalcs;

        mRefBases = refBases;
    }

    public MsiJitterCalcs msiJitterCalcs() { return mMsiJitterCalcs; }

    public static int calcEventPenalty(double numEvents, int readLength, double readMapQualEventsPenalty)
    {
        // penalise events as a function of read length
        return (int)round(max(0, numEvents - 1) / readLength * readMapQualEventsPenalty);
    }

    public static int modifiedMapQuality(
            final QualityConfig config, final BasePosition position, int mapQuality, int readLength, double readEvents, boolean isImproperPair)
    {
        if(isHighlyPolymorphic(position))
        {
            return min(HIGHLY_POLYMORPHIC_GENES_MAX_QUALITY, mapQuality - config.FixedMapQualPenalty);
        }

        int improperPairPenalty = isImproperPair ? config.ImproperPairPenalty : 0;
        int eventPenalty = calcEventPenalty(readEvents, readLength, config.ReadMapQualEventsPenalty);

        int modifiedMapQuality = mapQuality - config.FixedMapQualPenalty - improperPairPenalty - eventPenalty;

        return config.MapQualityRatioFactor > 0 ? min(MAX_MAP_QUALITY, modifiedMapQuality) : modifiedMapQuality;
    }

    public QualityScores calculateQualityScores(
            final ReadContextCounter readContextCounter, int readBaseIndex, final SAMRecord record, double numberOfEvents, double calcBaseQuality)
    {
        double baseQuality;

        boolean recalibrateBaseQuality;

        if(!isUltima())
        {
            recalibrateBaseQuality = !(readContextCounter.isIndel() || readContextCounter.artefactContext() != null);
        }
        else
        {
            recalibrateBaseQuality = readContextCounter.isSnv();
        }

        if(recalibrateBaseQuality)
        {
            if(isUltima())
            {
                double bqrQual = readContextCounter.qualCache().getQual(
                        UltimaUtils.maxRawQual(), SamRecordUtils.extractConsensusType(record), 0, !record.getReadNegativeStrandFlag());

                baseQuality = min(calcBaseQuality, bqrQual);
            }
            else
            {
                baseQuality = recalibratedBaseQuality(readContextCounter, readBaseIndex, record, readContextCounter.variant().ref().length());
            }
        }
        else
        {
            baseQuality = calcBaseQuality;
        }

        int mapQuality = record.getMappingQuality();
        boolean isImproperPair = isImproperPair(record);

        int modifiedMapQuality = modifiedMapQuality(
                mConfig, readContextCounter.variant(), mapQuality, record.getReadBases().length, numberOfEvents, isImproperPair);

        double modifiedBaseQuality = baseQuality;

        if(!readContextCounter.useMsiErrorRate())
            modifiedBaseQuality -= mConfig.BaseQualityFixedPenalty;

        int readEdgePenalty = readEdgeDistancePenalty(readContextCounter, readBaseIndex, record);
        modifiedBaseQuality -= readEdgePenalty;

        double modifiedQuality = max(0, min(modifiedMapQuality, modifiedBaseQuality));

        return new QualityScores(
                calcBaseQuality, baseQuality, max(0, modifiedMapQuality), max(0.0, modifiedBaseQuality), modifiedQuality);
    }

    public static boolean isHighBaseQual(final double baseQual)
    {
        return BaseQualAdjustment.isHighBaseQual((byte)baseQual, SEQUENCING_TYPE);
    }

    public static boolean isMediumBaseQual(final double baseQual)
    {
        return BaseQualAdjustment.isMediumBaseQual((byte)baseQual, SEQUENCING_TYPE);
    }

    public static boolean isImproperPair(final SAMRecord record) { return record.getReadPairedFlag() && !record.getProperPairFlag(); }

    public static double calculateBaseQuality(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        if(isUltima())
        {
            return readContextCounter.ultimaData().getQualModels().calculateQual(readContextCounter, readIndex, record);
        }

        byte artefactAdjustedQual = readContextCounter.artefactContext() != null ?
                readContextCounter.artefactContext().findApplicableBaseQual(record, readIndex) : INVALID_BASE_QUAL;

        if(readContextCounter.isIndel())
        {
            if(readContextCounter.qualCache().usesMsiIndelErrorQual() && artefactAdjustedQual != INVALID_BASE_QUAL)
            {
                // min the min of the two models
                return min(artefactAdjustedQual, readContextCounter.qualCache().msiIndelErrorQual());
            }

            double avgCoreQuality = isSbx() ?
                    SbxUtils.indelCoreQuality(readContextCounter.readContext(), record, readIndex) :
                    averageCoreQuality(readContextCounter.readContext(), record, readIndex);

            if(readContextCounter.qualCache().usesMsiIndelErrorQual())
                return min(avgCoreQuality, readContextCounter.qualCache().msiIndelErrorQual());
            else if(artefactAdjustedQual != INVALID_BASE_QUAL)
                return min(avgCoreQuality, artefactAdjustedQual);
            else
                return avgCoreQuality;
        }

        if(artefactAdjustedQual != INVALID_BASE_QUAL)
            return artefactAdjustedQual;

        if(readContextCounter.isSnv())
            return record.getBaseQualities()[readIndex];

        int varLength = readContextCounter.variant().ref().length();

        // take the minimum base qual for MNVs
        byte minBaseQual = record.getBaseQualities()[readIndex];

        for(int i = readIndex + 1; i < readIndex + varLength; ++i)
        {
            minBaseQual = minQual(minBaseQual, record.getBaseQualities()[i]);
        }

        return minBaseQual;
    }

    private double recalibratedBaseQuality(
            final ReadContextCounter readContextCounter, int startReadIndex, final SAMRecord record, int length)
    {
        ConsensusType consensusType = SamRecordUtils.extractConsensusType(record);
        boolean posOrientRead = !record.getReadNegativeStrandFlag();

        int dualBaseIndex = -1;
        if(isSbx() && consensusType == DUAL)
        {
            dualBaseIndex = SbxBamUtils.extractDuplexBaseIndex(record);

            if(!SbxBamUtils.inDuplexRegion(posOrientRead, dualBaseIndex, startReadIndex))
                consensusType = SINGLE;
        }

        if(readContextCounter.isSnv())
        {
            // simplified version of the MNV case below
            byte rawQuality = record.getBaseQualities()[startReadIndex];
            return readContextCounter.qualCache().getQual(rawQuality, consensusType, 0, posOrientRead);
        }

        // MNV case
        int maxIndex = min(startReadIndex + length, record.getBaseQualities().length) - 1;
        int maxLength = maxIndex - startReadIndex + 1;

        double quality = Integer.MAX_VALUE;
        for(int i = 0; i < maxLength; i++)
        {
            int readIndex = startReadIndex + i;
            byte rawQuality = record.getBaseQualities()[readIndex];

            double recalibratedQual = readContextCounter.qualCache().getQual(rawQuality, consensusType, i, posOrientRead);
            quality = min(quality, recalibratedQual);
        }

        return quality;
    }

    public byte[] getTrinucleotideContext(int refPosition)
    {
        return mRefBases.containsPosition(refPosition) ? mRefBases.trinucleotideContext(refPosition) : null;
    }

    public double lookupRecalibrateQuality(final byte[] trinucleotideContext, byte altBase, byte rawQuality, final ConsensusType readType)
    {
        if(rawQuality == 0)
            return 0; // never adjust a zero qual up

        if(mQualityRecalibrationMap == null)
            return rawQuality;

        return mQualityRecalibrationMap.getQualityAdjustment(trinucleotideContext[1], altBase, trinucleotideContext, rawQuality, readType);
    }

    public static double averageCoreQuality(final VariantReadContext readContext, final SAMRecord record, final int readVarIndex)
    {
        int readIndexStart = max(readVarIndex - readContext.leftCoreLength(), 0);
        int readIndexEnd = min(readVarIndex + readContext.rightCoreLength(), record.getReadBases().length - 1);

        int baseLength = readIndexEnd - readIndexStart + 1;

        if(baseLength <= 0)
            return 0;

        double quality = 0;

        for(int i = readIndexStart; i <= readIndexEnd; i++)
        {
            quality += record.getBaseQualities()[i];
        }

        return (int)round(quality / baseLength);
    }

    private int readEdgeDistancePenalty(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        int minDistance;

        if(readContextCounter.isIndel())
        {
            VariantReadContext readContext = readContextCounter.readContext();
            int lowerVarIndex = readIndex - (readContext.VarIndex - readContextCounter.matcher().altIndexLower());
            int upperVarIndex = readIndex + (readContextCounter.matcher().altIndexUpper() - readContext.VarIndex);

            minDistance = max(0, min(lowerVarIndex, record.getReadBases().length - 1 - upperVarIndex));
        }
        else
        {
            int upperReadIndex = readIndex + readContextCounter.variant().altLength() - 1; // for MNVs
            minDistance = max(0, min(readIndex, record.getReadBases().length - 1 - upperReadIndex));
        }

        if(minDistance <= 0)
            return READ_EDGE_PENALTY_0;
        else if(minDistance == 1)
            return READ_EDGE_PENALTY_1;
        else
            return 0;
    }
}
