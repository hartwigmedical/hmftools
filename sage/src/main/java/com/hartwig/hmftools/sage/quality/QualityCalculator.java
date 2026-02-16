package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.bam.ConsensusType.DUAL;
import static com.hartwig.hmftools.common.bam.ConsensusType.SINGLE;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.minQual;
import static com.hartwig.hmftools.sage.ReferenceData.isHighlyPolymorphic;
import static com.hartwig.hmftools.sage.SageCommon.isImproperPair;
import static com.hartwig.hmftools.sage.SageConfig.SEQUENCING_TYPE;
import static com.hartwig.hmftools.sage.SageConfig.isSbx;
import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.SageConstants.HIGHLY_POLYMORPHIC_GENES_MAX_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.MAX_MAP_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.READ_EDGE_PENALTY_0;
import static com.hartwig.hmftools.sage.SageConstants.READ_EDGE_PENALTY_1;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.JitterMatch;
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
        double nmProportion = (double) round(max(0, numEvents - 1)) / readLength;
        return (int) round(nmProportion * readMapQualEventsPenalty);
    }

    public QualityScores calculateQualityScores(
            final ReadContextCounter readContextCounter, int readVarIndex, final SAMRecord record, double numberOfEvents)
    {
        boolean usesMsiIndelErrorQual = readContextCounter.qualCache().usesMsiIndelErrorQual();
        boolean hasValidMsiJitterQuals = false;

        if(usesMsiIndelErrorQual)
        {
            hasValidMsiJitterQuals = JitterMatch.hasValidBaseQuals(
                    readContextCounter.readContext(), readContextCounter.matcher().altIndexUpper(), record, readVarIndex);
        }

        double seqTechBaseQuality = calcSequencingTechBaseQuality(readContextCounter, readVarIndex, record, hasValidMsiJitterQuals);

        if(seqTechBaseQuality == INVALID_BASE_QUAL)
            return QualityScores.INVALID_QUAL_SCORES;

        double recalibratedBaseQuality = seqTechBaseQuality;

        if(isUltima())
        {
            if(readContextCounter.isSnv())
            {
                double bqrQual = readContextCounter.qualCache().getQual(
                        UltimaUtils.maxRawQual(), SamRecordUtils.extractConsensusType(record), 0, !record.getReadNegativeStrandFlag());

                recalibratedBaseQuality = min(seqTechBaseQuality, bqrQual);
            }
        }
        else
        {
            if(!readContextCounter.isIndel() && readContextCounter.artefactContext() == null)
            {
                recalibratedBaseQuality = recalibratedBaseQuality(readContextCounter, readVarIndex, record);
                int varLength = readContextCounter.variant().ref().length();

                if(isSbx())
                {
                    if(readContextCounter.readContext().hasIndelInCore()
                    && seqTechBaseQuality < minBaseQualAcrossRange(readVarIndex, readVarIndex + varLength - 1, record))
                    {
                        // don't boost recalibrated qual of SBX artefacts
                        recalibratedBaseQuality = min(recalibratedBaseQuality, seqTechBaseQuality);
                    }
                }
            }
        }

        int mapQuality = record.getMappingQuality();
        boolean isImproperPair = isImproperPair(record);

        int penalisedMapQuality = calcMapQuality(
                mConfig, readContextCounter.variant(), mapQuality, record.getReadBases().length, numberOfEvents, isImproperPair);

        double penalisedBaseQuality = recalibratedBaseQuality;

        if(!readContextCounter.useMsiErrorRate())
            penalisedBaseQuality -= mConfig.BaseQualityFixedPenalty;

        int readEdgePenalty = readEdgeDistancePenalty(readContextCounter, readVarIndex, record);
        penalisedBaseQuality -= readEdgePenalty;

        double combinedQuality = max(0, min(penalisedMapQuality, penalisedBaseQuality));

        boolean isMediumQual = usesMsiIndelErrorQual ? !hasValidMsiJitterQuals : isMediumBaseQual(seqTechBaseQuality);

        return new QualityScores(
                seqTechBaseQuality, recalibratedBaseQuality, max(0, penalisedMapQuality), max(0.0, penalisedBaseQuality), combinedQuality,
                isMediumQual);
    }

    private static int calcMapQuality(
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

    public static boolean isHighBaseQual(final double baseQual)
    {
        return BaseQualAdjustment.isHighBaseQual((byte)baseQual, SEQUENCING_TYPE);
    }

    public static boolean isMediumBaseQual(final double baseQual)
    {
        return BaseQualAdjustment.isMediumBaseQual((byte)baseQual, SEQUENCING_TYPE);
    }

    private static double calcSequencingTechBaseQuality(
            final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record, boolean hasValidMsiJitterQuals)
    {
        if(isUltima())
        {
            return readContextCounter.ultimaData().getQualModels().calculateQual(readContextCounter, readIndex, record);
        }

        byte artefactAdjustedQual = readContextCounter.artefactContext() != null ?
                readContextCounter.artefactContext().findApplicableBaseQual(record, readIndex) : INVALID_BASE_QUAL;

        if(readContextCounter.isIndel())
        {
            double msiIndelErrorQualToUse = readContextCounter.qualCache().msiIndelErrorQual();

            if(readContextCounter.qualCache().usesMsiIndelErrorQual())
            {
                if(!hasValidMsiJitterQuals)
                {
                    msiIndelErrorQualToUse /= 2;
                }
                if(readContextCounter.allowUncertainCoreBases())
                    return msiIndelErrorQualToUse;
            }

            if(readContextCounter.qualCache().usesMsiIndelErrorQual() && artefactAdjustedQual != INVALID_BASE_QUAL)
            {
                // min the min of the two models
                return min(artefactAdjustedQual, msiIndelErrorQualToUse);
            }

            double avgCoreQuality = isSbx() ?
                    SbxUtils.indelCoreQuality(readContextCounter.readContext(), record, readIndex) :
                    averageCoreQuality(readContextCounter.readContext(), record, readIndex);

            if(readContextCounter.qualCache().usesMsiIndelErrorQual())
                return min(avgCoreQuality, msiIndelErrorQualToUse);
            else if(artefactAdjustedQual != INVALID_BASE_QUAL)
                return min(avgCoreQuality, artefactAdjustedQual);
            else
                return avgCoreQuality;
        }

        if(artefactAdjustedQual != INVALID_BASE_QUAL)
            return artefactAdjustedQual;

        if(isSbx() && readContextCounter.readContext().hasIndelInCore())
            return SbxUtils.indelCoreQuality(readContextCounter.readContext(), record, readIndex);

        int varLength = readContextCounter.variant().ref().length();
        return minBaseQualAcrossRange(readIndex, readIndex + varLength - 1, record);
    }

    public static double minBaseQualAcrossRange(final int startIndex, final int endIndex, final SAMRecord record)
    {
        byte minBaseQual = record.getBaseQualities()[startIndex];
        for(int i = startIndex + 1; i <= endIndex; ++i)
            minBaseQual = minQual(minBaseQual, record.getBaseQualities()[i]);
        return minBaseQual;
    }

    private double recalibratedBaseQuality(final ReadContextCounter readContextCounter, int startReadIndex, final SAMRecord record)
    {
        int variantLength = readContextCounter.variant().refLength();

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
        int maxIndex = min(startReadIndex + variantLength, record.getBaseQualities().length) - 1;
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
        if(rawQuality <= BaseQualAdjustment.BASE_QUAL_MINIMUM)
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
